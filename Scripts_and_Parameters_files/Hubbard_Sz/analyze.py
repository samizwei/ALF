#!/usr/bin/env python3
"""Analyse ALF runs of the Hubbard_Sz and TwoBand_Spin models.

Works on the text output of ALF after ana.out has run (run.sh does that).

Usage
-----
  python3 analyze.py RUN_DIR                 summary: scalars, sign, structure factors
  python3 analyze.py RUN_DIR --plot          + bins.png (thermalisation) and Sq.png (S(q) maps)
  python3 analyze.py RUN_DIR --all-q         print every q of the structure factors
  python3 analyze.py RUN1 RUN2 RUN3 --dtau   Dtau -> 0 extrapolation O(Dtau) = O_0 + b Dtau^2
                                              (runs that differ only in Dtau)
Only numpy is required; matplotlib only for --plot.

File formats used (ALF text output)
-----------------------------------
  X_scal    raw bins, one line per bin:  n  (Re,Im) (Re,Im) ...  sign
  X_scalJ   jackknife result: lines "i  mean  error" (last i = average sign)
  X_eqJK    one block per q:  "kx ky" then norb^2 lines "o1 o2 Re errRe Im errIm"
            (Fourier transform of the background-subtracted correlation;
             for SpinZ etc. with <S> = 0 this is S(q) = 1/N sum_ij e^{iq(ri-rj)} <S_i S_j>)
"""
import argparse
import os
import re
import sys

import numpy as np

SCALARS = ["Ener", "Kin", "Pot", "Part", "Part_cd"]
CORRS = ["SpinX", "SpinY", "SpinZ", "Pair"]          # one orbital per unit cell
LABELS = {"Ener": ["Ener"], "Kin": ["Kin"], "Pot": ["Pot"], "Part": ["Part"],
          "Part_cd": ["N_c", "N_d"]}


# ----------------------------------------------------------------------------- reading
def read_parameters(run):
    """Return the parameters used by this run as a dict (lower-case keys).

    Only the namelists that ALF reads for this model are used (VAR_ham_name, VAR_lattice,
    VAR_Model_Generic, VAR_<ham_name>, VAR_QMC, VAR_errors): a parameter file may contain
    e.g. Beta in several model namelists.
    """
    nml, cur = {}, None
    with open(os.path.join(run, "parameters")) as f:
        for line in f:
            line = line.split("!")[0].strip()
            if line.startswith("&"):
                cur = line[1:].split()[0].lower() if len(line) > 1 else None
                nml[cur] = {}
                continue
            if line.startswith("/"):
                cur = None
                continue
            m = re.match(r"(\w+)\s*=\s*(\S+)", line)
            if cur is None or not m:
                continue
            val = m.group(2).rstrip(",").strip('"\'')
            try:
                val = float(val.replace("d", "e").replace("D", "e"))
            except ValueError:
                pass
            nml[cur][m.group(1).lower()] = val
    ham = str(nml.get("var_ham_name", {}).get("ham_name", "")).lower()
    par = {}
    for name in ("var_ham_name", "var_lattice", "var_model_generic", "var_" + ham, "var_qmc", "var_errors"):
        par.update(nml.get(name, {}))
    return par


def read_scalJ(run, name):
    """Return (values, errors, sign, sign_err) from RUN/name_scalJ, or None."""
    path = os.path.join(run, name + "_scalJ")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 3 and p[0].isdigit():
                rows.append((float(p[1]), float(p[2])))
    rows = np.array(rows)
    return rows[:-1, 0], rows[:-1, 1], rows[-1, 0], rows[-1, 1]


def read_scal_bins(run, name):
    """Raw bins of RUN/name_scal: array (Nbins, Nobs) of Re(<O>), and the sign per bin."""
    path = os.path.join(run, name + "_scal")
    vals, sign = [], []
    with open(path) as f:
        for line in f:
            tup = re.findall(r"\(\s*([-+0-9.Ee]+)\s*,\s*([-+0-9.Ee]+)\s*\)", line)
            if not tup:
                continue
            vals.append([float(re_) for re_, _ in tup])
            sign.append(float(line.split(")")[-1]))
    vals, sign = np.array(vals), np.array(sign)
    return vals / sign[:, None], sign          # bins store <O s>/... -> divide by sign


def read_eqJK(run, name):
    """Return q (Nq,2), val (Nq,norb,norb), err (Nq,norb,norb) from RUN/name_eqJK."""
    path = os.path.join(run, name + "_eqJK")
    if not os.path.exists(path):
        return None
    qs, blocks = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                qs.append([float(p[0]), float(p[1])])
                blocks.append([])
            elif len(p) == 6:
                blocks[-1].append([int(p[0]), int(p[1]), float(p[2]), float(p[3])])
    norb = int(round(np.sqrt(len(blocks[0]))))
    val = np.zeros((len(qs), norb, norb))
    err = np.zeros_like(val)
    for iq, b in enumerate(blocks):
        for o1, o2, v, e in b:
            val[iq, o1 - 1, o2 - 1], err[iq, o1 - 1, o2 - 1] = v, e
    return np.array(qs), val, err


def q_index(qs, target):
    """Index of the q point equal to target modulo 2 pi."""
    d = np.mod(qs - np.array(target) + np.pi, 2 * np.pi) - np.pi
    i = np.argmin(np.linalg.norm(d, axis=1))
    return i if np.linalg.norm(d[i]) < 1e-6 else None


# ----------------------------------------------------------------------------- summary
def collect(run):
    """Dictionary label -> (mean, error) of the main observables of one run."""
    out = {}
    for name in SCALARS:
        r = read_scalJ(run, name)
        if r is None:
            continue
        v, e, s, se = r
        for lab, vi, ei in zip(LABELS[name], v, e):
            out[lab] = (vi, ei)
        if name == "Ener":
            out["Avg. sign"] = (s, se)
    for name in CORRS:
        r = read_eqJK(run, name)
        if r is None:
            continue
        qs, val, err = r
        for qlab, q in (("(pi,pi)", (np.pi, np.pi)), ("(0,0)", (0.0, 0.0))):
            i = q_index(qs, q)
            if i is not None:
                out["%s %s" % (name, qlab)] = (val[i, 0, 0], err[i, 0, 0])
    return out


def fmt(v, e):
    """Value with error in parentheses, e.g. -31.8098(72)."""
    if e <= 0 or not np.isfinite(e):
        return "%.6f" % v
    nd = max(0, 1 - int(np.floor(np.log10(e))))
    return "%.*f(%d)" % (nd, v, round(e * 10 ** nd))


def summary(run, all_q):
    par = read_parameters(run)
    print("Run: %s" % os.path.abspath(run))
    keys = ["ham_name", "ham_int", "l1", "l2", "ham_u", "beta", "dtau", "ham_chem"]
    print("  " + ", ".join("%s = %s" % (k, par[k]) for k in keys if k in par))
    r = read_scalJ(run, "Ener")
    if r is not None:
        with open(os.path.join(run, "Ener_scalJ")) as f:
            print("  " + f.readline().strip())
    print()
    print("  %-18s %s" % ("observable", "mean(error)"))
    for lab, (v, e) in collect(run).items():
        print("  %-18s %s" % (lab, fmt(v, e)))
    info = os.path.join(run, "info")
    if os.path.exists(info):
        for line in open(info):
            if "Precision Green" in line:
                print("\n  " + " ".join(line.split()))
    if all_q:
        for name in CORRS:
            r = read_eqJK(run, name)
            if r is None:
                continue
            qs, val, err = r
            print("\n  %s(q):      kx/pi    ky/pi   value(error)" % name)
            for q, v, e in zip(qs, val[:, 0, 0], err[:, 0, 0]):
                print("              %6.3f   %6.3f   %s" % (q[0] / np.pi, q[1] / np.pi, fmt(v, e)))


# ----------------------------------------------------------------------------- plots
def plots(run):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    par = read_parameters(run)
    n_skip = int(par.get("n_skip", 0))
    names = [n for n in SCALARS if os.path.exists(os.path.join(run, n + "_scal"))]
    fig, axes = plt.subplots(len(names), 1, figsize=(6, 1.8 * len(names)), sharex=True)
    for ax, name in zip(np.atleast_1d(axes), names):
        vals, _ = read_scal_bins(run, name)
        for k, lab in enumerate(LABELS[name]):
            ax.plot(np.arange(1, len(vals) + 1), vals[:, k], ".-", lw=0.8, label=lab)
        ax.axvspan(0.5, n_skip + 0.5, color="0.85", lw=0)
        ax.set_ylabel(name)
        if len(LABELS[name]) > 1:
            ax.legend(fontsize=8)
    np.atleast_1d(axes)[-1].set_xlabel("bin  (grey: discarded, n_skip = %d)" % n_skip)
    fig.tight_layout()
    fig.savefig(os.path.join(run, "bins.png"), dpi=150)
    print("  wrote %s" % os.path.join(run, "bins.png"))

    corrs = [n for n in CORRS if os.path.exists(os.path.join(run, n + "_eqJK"))]
    fig, axes = plt.subplots(1, len(corrs), figsize=(3.2 * len(corrs), 3))
    for ax, name in zip(np.atleast_1d(axes), corrs):
        qs, val, _ = read_eqJK(run, name)
        q = np.mod(qs + 1e-9, 2 * np.pi) - 1e-9          # map to [0, 2pi)
        kx, ky = np.unique(np.round(q[:, 0], 8)), np.unique(np.round(q[:, 1], 8))
        grid = np.full((len(ky), len(kx)), np.nan)
        for (a, b), v in zip(np.round(q, 8), val[:, 0, 0]):
            grid[np.searchsorted(ky, b), np.searchsorted(kx, a)] = v
        dx = kx[1] - kx[0] if len(kx) > 1 else 2 * np.pi
        dy = ky[1] - ky[0] if len(ky) > 1 else 2 * np.pi
        im = ax.imshow(grid, origin="lower", cmap="viridis",
                       extent=[(kx[0] - dx / 2) / np.pi, (kx[-1] + dx / 2) / np.pi,
                               (ky[0] - dy / 2) / np.pi, (ky[-1] + dy / 2) / np.pi])
        ax.set_title(name + "(q)")
        ax.set_xlabel(r"$q_x/\pi$")
        ax.set_ylabel(r"$q_y/\pi$")
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(run, "Sq.png"), dpi=150)
    print("  wrote %s" % os.path.join(run, "Sq.png"))


# ----------------------------------------------------------------------------- Dtau
def dtau_extrapolation(runs):
    dt = np.array([read_parameters(r)["dtau"] for r in runs])
    data = [collect(r) for r in runs]
    order = np.argsort(dt)[::-1]
    dt, data = dt[order], [data[i] for i in order]
    print("Weighted fit  O(Dtau) = O_0 + b Dtau^2  over Dtau = %s\n" % ", ".join("%g" % d for d in dt))
    print("  %-18s" % "observable" + "".join("%18s" % ("Dtau=%g" % d) for d in dt)
          + "%18s %8s" % ("O_0", "chi2/dof"))
    for lab in data[0]:
        if lab == "Avg. sign" or not all(lab in d for d in data):
            continue
        y = np.array([d[lab][0] for d in data])
        e = np.array([d[lab][1] for d in data])
        line = "  %-18s" % lab + "".join("%18s" % fmt(a, b) for a, b in zip(y, e))
        if len(dt) >= 2 and np.all(e > 0):
            A = np.vstack([np.ones_like(dt), dt ** 2]).T / e[:, None]
            coef, *_ = np.linalg.lstsq(A, y / e, rcond=None)
            cov = np.linalg.inv(A.T @ A)
            chi2 = np.sum(((A @ coef) - y / e) ** 2)
            dof = len(dt) - 2
            line += "%18s %8s" % (fmt(coef[0], np.sqrt(cov[0, 0])),
                                  "%.2f" % (chi2 / dof) if dof > 0 else "-")
        print(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run directories (after ana.out / run.sh)")
    ap.add_argument("--plot", action="store_true", help="write bins.png and Sq.png into each run directory")
    ap.add_argument("--all-q", action="store_true", help="print the structure factors at every q")
    ap.add_argument("--dtau", action="store_true", help="Dtau^2 extrapolation over the given runs")
    a = ap.parse_args()
    for r in a.runs:
        if not os.path.exists(os.path.join(r, "Ener_scalJ")):
            sys.exit("%s: no Ener_scalJ. Run ana.out first (run.sh does it)." % r)
    if a.dtau:
        dtau_extrapolation(a.runs)
        return
    for r in a.runs:
        summary(r, a.all_q)
        if a.plot:
            plots(r)
        print()


if __name__ == "__main__":
    main()
