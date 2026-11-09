#!/usr/bin/env bash
set -euo pipefail
code_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
exec "${RSCRIPT_BIN:-Rscript}" --vanilla "$code_dir/run.R" "$@"
