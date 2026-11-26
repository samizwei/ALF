# Hubbard_Sz: running the S^z S^z model

Two-band model of Christensen *et al.*, PRL **125**, 247001 (2020), Eq. (3):

    H = H_0 - U Σ_i (S^z_i)²,   S^z_i = Σ_ab c†_ia σ^z_ab d_ib + h.c.

Code: `Prog/Hamiltonians/Hamiltonian_Hubbard_Sz_smod.F90`.
Validation against ED: `Infos/Test_results_Hubbard_Sz.tex` (project folder).

Full manual (installation, data layout, analysis): `Infos_SZ/ALF_SZ_Manual.pdf`.
XY and Heisenberg models: `../TwoBand_Spin/` (same `run.sh`).

## Files

| file | purpose |
|---|---|
| `parameters` | template for your runs; lines marked `<-- SET` are the ones you usually change |
| `parameters_test_2x2` | 2×2 test with known exact (ED) results |
| `seeds` | random seeds (one line per MPI rank) |
| `run.sh` | sets up a run directory, runs ALF, runs the analysis, prints a summary |
| `out_to_in.sh` | used by `run.sh` to continue a run from its last configuration |
| `analyze.py` | summary table, plots (`--plot`) and Δτ² extrapolation (`--dtau`); see `Infos_SZ/ALF_SZ_Manual.pdf` |

## 0. Build ALF (once, or after changing Fortran code)

    cd ALF_SZ/ALF
    source configure.sh GNU noMPI     # or: GNU MPI  for parallel runs
    make

## 1. Check that everything works (a few seconds)

    cd ALF_SZ/ALF/Scripts_and_Parameters_files/Hubbard_Sz
    ./run.sh ~/qmc_runs/test_2x2 parameters_test_2x2

Expected (ED, β = 2, U = 1, Ham_chem = 0.3):

| quantity | ED | QMC (this test, Δτ = 0.05) |
|---|---|---|
| Ener | −31.80783 | −31.8098(72) |
| Part | 8.27457 | 8.2675(27) |
| N_c | 2.38350 | 2.3771(26) |
| N_d | 5.89106 | 5.8903(20) |
| Avg. sign | 1 | 1 |

N_c is about 2.5σ off at Δτ = 0.05. That is Trotter error, and it goes away at smaller Δτ.

## 2. Your own runs

1. Copy the template, e.g. `cp parameters parameters_L8_U2`, and edit the `<-- SET` lines:
   L1, L2, Ham_U, Beta, Dtau, NSweep, NBin, n_skip.
   δ = 0.4 and μ = −2 are already set to the paper's values.
2. Start it in its own directory:

       ./run.sh ~/qmc_runs/L8_U2 parameters_L8_U2

3. To collect more statistics, run the same command again. The run continues from the
   last configuration and appends new bins.

With an MPI build: `NP=4 ./run.sh <dir> <parameters>`.

## 3. Output (inside the run directory)

| file | content |
|---|---|
| `Ener_scalJ`, `Kin_scalJ`, `Pot_scalJ`, `Part_scalJ` | mean and jackknife error; the second entry in `Ener_scalJ` is the average sign |
| `Part_cd_scalJ` | N_c (entry 1) and N_d (entry 2) |
| `SpinZ_eqJK`, `Pair_eqJK`, `Den_eqJK`, `Green_eqJK` | equal-time correlations in momentum space; one block per q = (kx, ky), value and error |
| `*_eqJR` | the same in real space |
| `*_tauJK` | time-displaced correlations (only with `Ltau = 1`) |
| `info` | run settings and numerical stability ("Precision Green"; should be ≪ 1e-6) |
| `run.log`, `ana.log` | program output |

`SpinZ_eqJK` at q = (π,π) is the S^z structure factor. `Pair_eqJK` at q = 0 is the pair
structure factor with P_i = 2(c↑c↓ − d↑d↓).

## Practical notes

- **Thermalization:** check that the first bins are not drifting, e.g.
  `grep -v '^#' Ener_scal | head`. If they are, increase `n_skip` and rerun only the analysis:
  `<ALF>/Analysis/ana.out *_scal *_eq`.
- **Δτ:** the Trotter error scales as Δτ². For final numbers, compare two Δτ values or
  extrapolate.
- **Cost:** it grows roughly as β·N³ with N = L1·L2. The 2×2 test takes seconds. The default
  template (6×6, β = 5, Δτ = 0.05, `Ltau = 0`) takes about 20 s per bin of 200 sweeps on
  this laptop; `Ltau = 1` makes it noticeably slower.
- `run.sh` always uses the ALF folder two levels above it. A global `ALF_DIR` variable,
  which on this machine points to a different ALF installation, is ignored.
