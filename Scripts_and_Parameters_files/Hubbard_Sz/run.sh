#!/bin/bash
# Run the Hubbard_Sz (Ising, S^z S^z) model with ALF and analyse the output.
#
# Usage:
#   ./run.sh <run_dir> [parameters_file]
#
#   <run_dir>          directory for this run (created if it does not exist)
#   [parameters_file]  default: ./parameters next to this script
#
# If <run_dir> already contains a finished run (confout_*), the run is continued:
# the last configuration is used as the start and new bins are appended.
#
# The ALF executable is always taken from the ALF folder that contains this script
# (two levels up), not from a global ALF_DIR variable.
# Environment variable:
#   NP        number of MPI processes (only if ALF was built with MPI; default: serial)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALF_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXE="$ALF_DIR/Prog/ALF.out"
ANA="$ALF_DIR/Analysis/ana.out"

if [ $# -lt 1 ]; then
  sed -n '2,16p' "$0"
  exit 1
fi
RUN_DIR="$1"
PARAMS="${2:-$SCRIPT_DIR/parameters}"

for f in "$EXE" "$ANA" "$PARAMS"; do
  if [ ! -e "$f" ]; then
    echo "Not found: $f"
    [ "$f" = "$EXE" ] && echo "Build ALF first:  cd $ALF_DIR && source configure.sh GNU noMPI && make"
    exit 1
  fi
done

echo "Using ALF in $ALF_DIR"
mkdir -p "$RUN_DIR"
PARAMS="$(cd "$(dirname "$PARAMS")" && pwd)/$(basename "$PARAMS")"
cd "$RUN_DIR"

if ls confout_* >/dev/null 2>&1; then
  echo "Continuing previous run in $(pwd)"
  bash "$SCRIPT_DIR/out_to_in.sh"
  if ! cmp -s "$PARAMS" parameters; then
    echo "Note: parameters differ from the previous run; using the new file."
  fi
else
  echo "Starting new run in $(pwd)"
fi
cp "$PARAMS" parameters
[ -e seeds ] || cp "$SCRIPT_DIR/seeds" seeds

echo "Running ALF ... (log: $(pwd)/run.log)"
START=$(date +%s)
if [ -n "$NP" ]; then
  mpirun -np "$NP" "$EXE" > run.log 2>&1
else
  "$EXE" > run.log 2>&1
fi
echo "ALF finished in $(( $(date +%s) - START )) s"

echo "Analysing ..."
EQ=$(ls *_eq 2>/dev/null || true)
TAU=$(ls *_tau 2>/dev/null || true)
"$ANA" *_scal $EQ $TAU > ana.log 2>&1

# Short summary: scalar observables (mean, jackknife error) and average sign
echo
echo "Observable           mean              error"
for f in *_scalJ; do
  name="${f%_scalJ}"
  awk -v n="$name" '!/^#/ && NF==3 && $1>=1 {
        if (n=="Part_cd") lab = ($1==1 ? "N_c" : ($1==2 ? "N_d" : ""));
        else              lab = ($1==1 ? n : "");
        if (lab != "") printf "%-12s %16.8f %16.8f\n", lab, $2, $3;
        else if (n=="Ener" && !s) {printf "%-12s %16.8f %16.8f\n", "Avg. sign", $2, $3; s=1}
      }' "$f"
done
grep -E "Precision" info | sed 's/^ */  /'
echo
echo "Correlations in momentum space: *_eqJK (and *_tauJK if Ltau = 1)"
