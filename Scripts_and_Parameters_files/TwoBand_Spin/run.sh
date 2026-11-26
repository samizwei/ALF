#!/bin/bash
# Run the TwoBand_Spin model (Ising / XY / Heisenberg) with ALF and analyse the output.
#
# Usage:
#   ./run.sh <run_dir> [parameters_file]
#
#   [parameters_file]  default: ./parameters next to this script
#
# This is a thin wrapper around ../Hubbard_Sz/run.sh, which is model independent
# (the model is chosen by ham_name in the parameters file). See that script for details.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ $# -lt 1 ]; then sed -n '2,10p' "$0"; exit 1; fi
exec "$SCRIPT_DIR/../Hubbard_Sz/run.sh" "$1" "${2:-$SCRIPT_DIR/parameters}"
