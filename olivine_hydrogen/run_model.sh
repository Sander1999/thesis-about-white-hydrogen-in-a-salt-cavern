#!/usr/bin/env bash
set -euo pipefail
MODEL_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODEL_PYTHON="${OLIVINE_PYTHON:-/Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python}"
MODEL_SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/olivine-run.XXXXXXXX")"
trap 'rm -rf -- "$MODEL_SCRATCH"' EXIT
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$MODEL_SCRATCH/matplotlib"
export IPYTHONDIR="$MODEL_SCRATCH/ipython"
export JUPYTER_RUNTIME_DIR="$MODEL_SCRATCH/jupyter"
export OMP_NUM_THREADS=2
cd "$MODEL_DIR/.."
case "${1:-all}" in
  data) "$MODEL_PYTHON" -B -m olivine_hydrogen.data.process_xrf ;;
  batch) "$MODEL_PYTHON" -B -m olivine_hydrogen.run_batch ;;
  feedback) "$MODEL_PYTHON" -B -m olivine_hydrogen.feedback.run ;;
  darts) "$MODEL_PYTHON" -B -m olivine_hydrogen.darts.run --suite ;;
  check) "$MODEL_PYTHON" -B -m unittest discover -s olivine_hydrogen/tests -v
         "$MODEL_PYTHON" -B -m unittest olivine_hydrogen.feedback.test_model -v
         "$MODEL_PYTHON" -B -m olivine_hydrogen.darts.test_model ;;
  notebook) "$MODEL_PYTHON" -B -m olivine_hydrogen.execute_notebook ;;
  all) "$MODEL_PYTHON" -B -m olivine_hydrogen.data.process_xrf
       "$MODEL_PYTHON" -B -m olivine_hydrogen.run_batch
       "$MODEL_PYTHON" -B -m olivine_hydrogen.feedback.run
       "$MODEL_PYTHON" -B -m olivine_hydrogen.darts.run --suite ;;
  *) echo 'Usage: run_model.sh [all|data|batch|feedback|darts|check|notebook]' >&2; exit 2 ;;
esac
