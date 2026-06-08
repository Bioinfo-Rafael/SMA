#!/usr/bin/env bash
# End-to-end pipeline. Override any hyperparameter via env vars, e.g.:
#     MAX_EPOCHS=200 RESOLUTIONS=0.2,0.5,1.0 ./run.sh
#     ./run.sh 03 04 05            # re-run only specific steps
set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-$PWD/.venv/bin/python}"

ALL=(01_concat 02_qc 03_scvi 04_umap 05_cluster 06_annotate 07_plot)

if [[ $# -eq 0 ]]; then
    STEPS=("${ALL[@]}")
else
    STEPS=()
    for arg in "$@"; do
        match=""
        for s in "${ALL[@]}"; do
            if [[ "$s" == "$arg"* ]]; then match="$s"; break; fi
        done
        [[ -z "$match" ]] && { echo "unknown step: $arg"; exit 2; }
        STEPS+=("$match")
    done
fi

for s in "${STEPS[@]}"; do
    echo
    echo "=================  $s  ================="
    "$PY" "scripts/${s}.py"
done
echo
echo "done."
