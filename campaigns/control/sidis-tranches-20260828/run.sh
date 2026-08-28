#!/usr/bin/env bash
set -euo pipefail

if [[ -r /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
elif [[ -r /opt/homebrew/opt/modules/init/bash ]]; then
  source /opt/homebrew/opt/modules/init/bash
else
  echo "error: no supported Environment Modules initialization was found" >&2
  exit 2
fi
module purge
module load herwig/pol

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CONTROL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$CONTROL_DIR/controller.py" "$@"
