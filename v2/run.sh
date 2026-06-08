#!/usr/bin/env bash
# v2 preprocessing driver. Runs steps 00..05 in order.
# Step 06 is a template and is NOT run automatically.
#
#   ./run.sh                 # all steps 00-05
#   ./run.sh 01 02 03        # only these steps
#   PYTHON=../v1/.venv/bin/python ./run.sh   # reuse the v1 venv
#
# Override the interpreter with the PYTHON env var. By default it looks for a
# local .venv, then falls back to `python3`.
set -euo pipefail
cd "$(dirname "$0")"

if [[ -n "${PYTHON:-}" ]]; then
    PY="$PYTHON"
elif [[ -x ".venv/bin/python" ]]; then
    PY=".venv/bin/python"
else
    PY="python3"
fi

ALL=(00_make_manifest 01_download_geo_supplement 02_extract_archives \
     03_load_to_anndata 04_inspect_anndata_columns 05_curate_and_save_h5ad)

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

echo "interpreter: $PY"
for s in "${STEPS[@]}"; do
    echo
    echo "=================  $s  ================="
    "$PY" "scripts/${s}.py"
done
echo
echo "done. (step 06_read_saved_h5ad_template.py is a template; run manually)"
