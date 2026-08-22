#!/usr/bin/env bash
set -euo pipefail

source /etc/profile.d/modules.sh
module purge
module load herwig/pol

THEPEG_OVERLAY=/home/apapaefs/Projects/Herwig/production/star510-bloch-ad0c5486/prefix
export LD_LIBRARY_PATH="$THEPEG_OVERLAY/lib/ThePEG:${LD_LIBRARY_PATH:-}"
export THEPEG_LIBRARY_PATH="$THEPEG_OVERLAY/lib/ThePEG:${THEPEG_LIBRARY_PATH:-}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CONTROL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$CONTROL_DIR/controller.py" "$@"
