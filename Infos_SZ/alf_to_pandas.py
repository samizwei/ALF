#!/usr/bin/env python3
"""Load the text output of ALF run directories into pandas DataFrames.

Works for every ALF Hamiltonian (Hubbard_Sz, TwoBand_Spin, the standard Hubbard
model, ...), for ALF builds without HDF5 (text output). The structure of the
tables is explained in ALF_pandas_guide.pdf in this folder.

Quick start (Python)
--------------------
    import sys; sys.path.insert(0, "/path/to/ALF_SZ/ALF/Infos_SZ")
    import alf_to_pandas as ap

    d = ap.load_run("/path/to/run_dir")          # one run
    d.scal                                        # ana.out results of scalar observables
    d.eq.query("obs == 'SpinZ'")                  # S(q) with error bars

    D = ap.load_runs(["run_U1", "run_U2"])        # several runs, stacked, column "run"
    D.scal.merge(D.params[["run", "ham_u"]], on="run")

Command line
------------
    python3 alf_to_pandas.py RUN_DIR [RUN_DIR ...]               # print an overview
    python3 alf_to_pandas.py RUN_DIR ... --pickle all.pkl        # save (pandas.read_pickle)
    python3 alf_to_pandas.py RUN_DIR ... --csv out_dir           # one CSV file per table

Requirements: numpy, pandas. ana.out must have been run in each run directory
(run.sh does that) for the analysed tables (scal, eq, eq_r, tau, chi); the raw-bin
tables only need the output of ALF.out.
"""
import argparse
import os
import re
from dataclasses import dataclass, field, fields

import numpy as np
import pandas as pd

__all__ = ["load_run", "load_runs", "jackknife", "jackknife_bins",
           "read_parameters", "ALFData"]

# Names of the components of vector-valued scalar observables
COMPONENT_LABELS = {"Part_cd": ["N_c", "N_d"]}

# Names of the orbitals of the unit cell for the multi-orbital observables (Green, Den)
ORBITAL_LABELS = {
    "hubbard_sz": ["c", "d"],
    "twoband_spin": ["c_up", "c_dn", "d_up", "d_dn"],
}

_FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?|NaN|nan|Infinity|-Infinity"
_CPLX = re.compile(r"\(\s*(" + _FLOAT + r")\s*,\s*(" + _FLOAT + r")\s*\)")


# ============================================================================ container
@dataclass
class ALFData:
    """All tables of one run (load_run) or of several runs (load_runs).

    Every table has a column "run" (the run directory name), so tables of
    different runs can be concatenated and merged with `params`.
    """
    params: pd.DataFrame = field(default_factory=pd.DataFrame)     # one row per run
    info: pd.DataFrame = field(default_factory=pd.DataFrame)       # one row per run
    scal: pd.DataFrame = field(default_factory=pd.DataFrame)       # *_scalJ
    scal_bins: pd.DataFrame = field(default_factory=pd.DataFrame)  # *_scal
    eq: pd.DataFrame = field(default_factory=pd.DataFrame)         # *_eqJK
    eq_r: pd.DataFrame = field(default_factory=pd.DataFrame)       # *_eqJR
    eq_bins: pd.DataFrame = field(default_factory=pd.DataFrame)    # *_eq  (raw=True)
    eq_back: pd.DataFrame = field(default_factory=pd.DataFrame)    # *_eq  (raw=True)
    tau: pd.DataFrame = field(default_factory=pd.DataFrame)        # *_<kx>_<ky>/g_dat, *_R0/g_dat
    chi: pd.DataFrame = field(default_factory=pd.DataFrame)        # *_tauJK
    tau_bins: pd.DataFrame = field(default_factory=pd.DataFrame)   # *_tau (raw_tau=True)

    def tables(self):
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def __repr__(self):
        lines = ["ALFData with tables (rows x columns):"]
        for name, df in self.tables().items():
            if len(df):
                lines.append(f"  {name:10s} {df.shape[0]:8d} x {df.shape[1]:<3d}  {', '.join(map(str, df.columns))}")
            else:
                lines.append(f"  {name:10s}        -")
        return "\n".join(lines)

    def to_pickle(self, path):
        pd.to_pickle(self.tables(), path)

    def to_csv(self, directory):
        os.makedirs(directory, exist_ok=True)
        for name, df in self.tables().items():
            if len(df):
                df.to_csv(os.path.join(directory, name + ".csv"), index=False)


# ============================================================================ parameters
def _parse_value(s):
    s = s.strip().rstrip(",").strip()
    if len(s) >= 2 and s[0] in "\"'" and s[-1] == s[0]:
        return s[1:-1]
    low = s.lower()
    if low in (".t.", ".true.", "t", "true"):
        return True
    if low in (".f.", ".false.", "f", "false"):
        return False
    try:
        if re.fullmatch(r"[-+]?\d+", s):
            return int(s)
        return float(s.replace("d", "e").replace("D", "e"))
    except ValueError:
        return s


def read_namelists(run):
    """All namelists of RUN/parameters as {namelist_name: {key: value}} (lower-case)."""
    nml, current = {}, None
    with open(os.path.join(run, "parameters")) as f:
        for raw in f:
            line = raw.split("!")[0].strip()
            if not line:
                continue
            if line.startswith("&"):
                current = line[1:].split()[0].lower()
                nml[current] = {}
                continue
            if line.startswith("/"):
                current = None
                continue
            if current is None:
                continue
            for m in re.finditer(r"(\w+)\s*=\s*(\"[^\"]*\"|'[^']*'|[^,\s]+)", line):
                nml[current][m.group(1).lower()] = _parse_value(m.group(2))
    return nml


def read_parameters(run):
    """Parameters that this run actually used, as one flat dict with lower-case keys.

    Only the namelists read by ALF for this model are used: VAR_ham_name, VAR_lattice,
    VAR_Model_Generic, VAR_<ham_name>, VAR_QMC, VAR_errors. (A parameter file can
    contain e.g. "Beta" in several model namelists; only the right one is kept.)
    """
    nml = read_namelists(run)
    ham = str(nml.get("var_ham_name", {}).get("ham_name", "")).lower()
    out = {"ham_name": nml.get("var_ham_name", {}).get("ham_name", "")}
    for name in ("var_lattice", "var_model_generic", "var_" + ham, "var_qmc", "var_errors"):
        out.update(nml.get(name, {}))
    return out


def read_info(run):
    """Main entries of RUN/info for the last ALF call (a continued run appends to info)."""
    path = os.path.join(run, "info")
    if not os.path.exists(path):
        return {}
    text = open(path).read()
    segments = text.split("START TIME")
    last = segments[-1]
    out = {"n_alf_calls": max(1, len(segments) - 1)}
    pats = {
        "precision_green_mean": r"Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)",
        "precision_phase_max": r"Precision Phase, Max\s*:\s*(\S+)",
        "precision_tau_mean": r"Precision tau\s+Mean, Max\s*:\s*(\S+)\s+(\S+)",
        "acceptance": r"^\s*Acceptance\s*:\s*(\S+)",
        "cpu_time_s": r"CPU Time\s*:\s*(\S+)",
    }
    for key, pat in pats.items():
        m = re.search(pat, last, re.M)
        if m:
            out[key] = float(m.group(1))
            if key == "precision_green_mean":
                out["precision_green_max"] = float(m.group(2))
            if key == "precision_tau_mean":
                out["precision_tau_max"] = float(m.group(2))
    return out


# ============================================================================ scalars
def _labels(obs, n):
    if obs in COMPONENT_LABELS and len(COMPONENT_LABELS[obs]) == n:
        return COMPONENT_LABELS[obs]
    return [obs] if n == 1 else [f"{obs}_{i + 1}" for i in range(n)]


def _read_scal_bins(path, obs):
    rows = []
    with open(path) as f:
        for ib, line in enumerate(l for l in f if "(" in l):
            vals = [complex(float(a), float(b)) for a, b in _CPLX.findall(line)]
            sign = float(line.rsplit(")", 1)[1])
            for k, (lab, v) in enumerate(zip(_labels(obs, len(vals)), vals)):
                rows.append((obs, k + 1, lab, ib + 1, v.real, v.imag, sign))
    df = pd.DataFrame(rows, columns=["obs", "comp", "label", "bin", "Os_re", "Os_im", "sign"])
    df["value"] = df["Os_re"] / df["sign"]
    return df


def _read_scalJ(path, obs):
    nbins_eff = nbins = None
    rows = []
    with open(path) as f:
        for line in f:
            if "Effective" in line:
                nbins_eff, nbins = map(int, line.split()[-2:])
                continue
            p = line.split()
            if len(p) == 3 and p[0].lstrip("-").isdigit():
                rows.append((int(p[0]), float(p[1]), float(p[2])))
    *obs_rows, sign_row = rows
    labs = _labels(obs, len(obs_rows))
    df = pd.DataFrame([(obs, i, lab, m, e) for (i, m, e), lab in zip(obs_rows, labs)],
                      columns=["obs", "comp", "label", "mean", "error"])
    df["sign"], df["sign_error"] = sign_row[1], sign_row[2]
    df["n_bins"], df["n_bins_eff"] = nbins, nbins_eff
    return df


# ============================================================================ lattice observables
def _read_corr_JK(path, obs, coord=("kx", "ky")):
    """*_eqJK / *_eqJR: blocks 'x y' followed by 'o1 o2 Re errRe Im errIm'."""
    rows, x = [], None
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                x = (float(p[0]), float(p[1]))
            elif len(p) == 6:
                rows.append((obs, *x, int(p[0]), int(p[1]), *map(float, p[2:])))
    return pd.DataFrame(rows, columns=["obs", *coord, "o1", "o2", "re", "re_error", "im", "im_error"])


def _read_latt_bins(path, obs):
    """Raw bins of an equal-time (*_eq) or time-displaced (*_tau) lattice observable.

    Per bin: header 'sign Norb Nunit [Ntau dtau]', Norb background values <O_o>,
    then for every k: 'kx ky' and Ntau*Norb*Norb complex values (order tau, o1, o2).
    Returns (bins, back) DataFrames; values are <O s> (divide by sign for <O>).
    """
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    rows, back, i, ib = [], [], 0, 0
    while i < len(lines):
        h = lines[i].split()
        sign, norb, nunit = float(h[0]), int(h[1]), int(h[2])
        ntau = int(h[3]) if len(h) >= 5 else 1
        dtau = float(h[4]) if len(h) >= 5 else 0.0
        ib += 1
        i += 1
        for o in range(norb):
            a, b = _CPLX.findall(lines[i])[0]
            back.append((obs, ib, o + 1, float(a), float(b), sign))
            i += 1
        for _ in range(nunit):
            kx, ky = map(float, lines[i].split())
            i += 1
            vals = np.array([[float(a), float(b)] for l in lines[i:i + ntau * norb * norb]
                             for a, b in _CPLX.findall(l)]).reshape(ntau, norb, norb, 2)
            i += ntau * norb * norb
            for nt in range(ntau):
                for o1 in range(norb):
                    for o2 in range(norb):
                        rows.append((obs, ib, kx, ky, nt * dtau, o1 + 1, o2 + 1,
                                     vals[nt, o1, o2, 0], vals[nt, o1, o2, 1], sign))
    bins = pd.DataFrame(rows, columns=["obs", "bin", "kx", "ky", "tau", "o1", "o2", "Os_re", "Os_im", "sign"])
    if ntau == 1:
        bins = bins.drop(columns="tau")
    backdf = pd.DataFrame(back, columns=["obs", "bin", "orb", "Os_re", "Os_im", "sign"])
    return bins, backdf


def _read_g_dat(path, obs, kx, ky):
    """ana.out output X_<kx>_<ky>/g_dat or X_R0/g_dat: header, then 'tau value error'."""
    with open(path) as f:
        h = f.readline().split()
        nt = int(h[0])
        channel = h[4] if len(h) > 4 else ""
        data = np.array([list(map(float, f.readline().split()[:3])) for _ in range(nt)])
    df = pd.DataFrame({"obs": obs, "kx": kx, "ky": ky, "tau": data[:, 0],
                       "value": data[:, 1], "error": data[:, 2]})
    df["channel"] = channel
    return df


def _read_tauJK(path, obs):
    d = np.loadtxt(path, ndmin=2)
    return pd.DataFrame({"obs": obs, "kx": d[:, 0], "ky": d[:, 1], "re": d[:, 2],
                         "re_error": d[:, 3], "im": d[:, 4], "im_error": d[:, 5]})


def _add_orbital_names(df, ham_name):
    labels = ORBITAL_LABELS.get(str(ham_name).lower())
    if not len(df) or "o1" not in df or labels is None:
        return df
    norb = df.groupby("obs")["o1"].transform("max")
    ok = norb == len(labels)
    for c in ("o1", "o2"):
        df[c + "_name"] = None
        df.loc[ok, c + "_name"] = [labels[o - 1] for o in df.loc[ok, c]]
        df.loc[~ok, c + "_name"] = df.loc[~ok, c].astype(str)
    return df


# ============================================================================ public loaders
def load_run(run, raw=True, raw_tau=False, name=None):
    """Read one ALF run directory.

    raw      read the raw bins of scalar and equal-time observables (scal_bins, eq_bins, eq_back)
    raw_tau  also read the raw bins of time-displaced observables (tau_bins; can be large)
    name     value of the "run" column (default: name of the directory)
    """
    run = os.path.abspath(run)
    name = name or os.path.basename(run.rstrip("/"))
    files = sorted(os.listdir(run))
    par = read_parameters(run)
    d = ALFData()
    d.params = pd.DataFrame([{"run": name, "path": run, **par}])
    d.info = pd.DataFrame([{"run": name, **read_info(run)}])

    def collect(suffix, reader, **kw):
        parts = [reader(os.path.join(run, f), f[: -len(suffix)], **kw)
                 for f in files if f.endswith(suffix) and os.path.isfile(os.path.join(run, f))]
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    d.scal = collect("_scalJ", _read_scalJ)
    d.eq = collect("_eqJK", _read_corr_JK)
    d.eq_r = collect("_eqJR", _read_corr_JK, coord=("rx", "ry"))
    d.chi = collect("_tauJK", _read_tauJK)
    if raw:
        d.scal_bins = collect("_scal", _read_scal_bins)
        parts = [_read_latt_bins(os.path.join(run, f), f[:-3]) for f in files if f.endswith("_eq")]
        if parts:
            d.eq_bins = pd.concat([p[0] for p in parts], ignore_index=True)
            d.eq_back = pd.concat([p[1] for p in parts], ignore_index=True)
    if raw_tau:
        parts = [_read_latt_bins(os.path.join(run, f), f[:-4])[0] for f in files if f.endswith("_tau")]
        if parts:
            d.tau_bins = pd.concat(parts, ignore_index=True)

    # time-displaced results of ana.out: directories <obs>_<kx>_<ky>/g_dat and <obs>_R0/g_dat
    tau_obs = {f[:-4] for f in files if f.endswith("_tau")}
    parts = []
    for f in files:
        gd = os.path.join(run, f, "g_dat")
        if not os.path.isfile(gd):
            continue
        m = re.fullmatch(r"(.+?)_(e_)?(\d+\.\d+)_(\d+\.\d+)", f)
        if f.endswith("_R0") and f[:-3] in tau_obs:
            parts.append(_read_g_dat(gd, f[:-3], np.nan, np.nan).assign(point="R0"))
        elif m and m.group(1) in tau_obs:
            parts.append(_read_g_dat(gd, m.group(1), float(m.group(3)), float(m.group(4)))
                         .assign(point="k" if not m.group(2) else "k_extended"))
    d.tau = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    for key, df in d.tables().items():
        if len(df) and key not in ("params", "info"):
            df.insert(0, "run", name)
            _add_orbital_names(df, par.get("ham_name", ""))
    return d


def load_runs(runs, **kw):
    """Read several run directories and stack all tables (column "run" tells them apart).

    `runs` is a list of directories or a glob pattern such as "~/QMC_runs/Heis_*".
    """
    if isinstance(runs, str):
        import glob
        runs = sorted(p for p in glob.glob(os.path.expanduser(runs)) if os.path.isdir(p))
    data = [load_run(r, **kw) for r in runs]
    out = ALFData()
    for key in out.tables():
        parts = [getattr(d, key) for d in data if len(getattr(d, key))]
        setattr(out, key, pd.concat(parts, ignore_index=True) if parts else pd.DataFrame())
    return out


# ============================================================================ jackknife
def jackknife(Os, sign, func=None, n_skip=0, n_rebin=1):
    """Jackknife mean and error, as in ana.out, from per-bin <O s> and <s>.

    Os     array (n_bins,) or (n_bins, m): sign-weighted bin averages (columns Os_re)
    sign   array (n_bins,): average sign per bin
    func   optional function of the m bin-averaged observables -> derived quantity,
           e.g. func=lambda O: O[0] / O[1]. Default: the observables themselves.
    n_skip bins discarded at the start; n_rebin consecutive bins combined into one.
    Returns (mean, error) as floats or arrays. The error follows the convention of
    ana.out (error^2 = N * variance of the jackknife estimates), so the numbers agree
    with the *_scalJ / *_eqJK files.
    """
    Os = np.asarray(Os, dtype=float)
    Os = Os[:, None] if Os.ndim == 1 else Os
    s = np.asarray(sign, dtype=float)[n_skip:]
    Os = Os[n_skip:]
    nb = len(s) // n_rebin
    if nb < 2:
        raise ValueError("need at least 2 (rebinned) bins after n_skip")
    Os = Os[: nb * n_rebin].reshape(nb, n_rebin, -1).sum(1)
    s = s[: nb * n_rebin].reshape(nb, n_rebin).sum(1)
    f = func or (lambda O: O)
    tot_O, tot_s = Os.sum(0), s.sum()
    jk = np.array([np.atleast_1d(f((tot_O - Os[j]) / (tot_s - s[j]))) for j in range(nb)])
    mean = jk.mean(0)
    # ALF convention (errors_mod.F90, ERRCALC_JS): err^2 = N * var_jk.  The textbook
    # jackknife uses (N - 1) instead of N; the difference is sqrt(N/(N-1)).
    err = np.sqrt(nb * ((jk - mean) ** 2).mean(0))
    if mean.size == 1:
        return float(mean[0]), float(err[0])
    return mean, err


def jackknife_bins(bins, n_skip=0, n_rebin=1, by=None):
    """Jackknife every group of a raw-bin table (scal_bins, eq_bins, tau_bins).

    Groups by all identifying columns (run, obs, comp/k/orbitals/tau) and returns a
    table with columns mean, error. Note: no background subtraction, i.e. for
    eq_bins this agrees with ana.out (table eq) at all q except q = 0 when <O> != 0.
    """
    if by is None:
        by = [c for c in bins.columns if c not in ("bin", "Os_re", "Os_im", "sign", "value")]
    out = []
    for key, g in bins.groupby(by, sort=False):
        g = g.sort_values("bin")
        m, e = jackknife(g["Os_re"].to_numpy(), g["sign"].to_numpy(), n_skip=n_skip, n_rebin=n_rebin)
        out.append((*key, m, e) if isinstance(key, tuple) else (key, m, e))
    return pd.DataFrame(out, columns=[*by, "mean", "error"])


# ============================================================================ CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="ALF run directories")
    ap.add_argument("--raw-tau", action="store_true", help="also read raw bins of time-displaced observables")
    ap.add_argument("--pickle", metavar="FILE", help="save all tables (dict of DataFrames) with pandas.to_pickle")
    ap.add_argument("--csv", metavar="DIR", help="save every table as DIR/<table>.csv")
    a = ap.parse_args()
    d = load_runs(a.runs, raw_tau=a.raw_tau)
    print(d)
    if len(d.scal):
        print("\nScalar observables (ana.out):")
        with pd.option_context("display.width", 160, "display.max_rows", 200):
            print(d.scal[["run", "label", "mean", "error", "sign", "n_bins_eff"]].to_string(index=False))
    if a.pickle:
        d.to_pickle(a.pickle)
        print(f"\nsaved: {a.pickle}   (load with: tables = pandas.read_pickle('{a.pickle}'))")
    if a.csv:
        d.to_csv(a.csv)
        print(f"\nsaved CSV files in {a.csv}")


if __name__ == "__main__":
    main()
