#!/usr/bin/env bash
# v2 DOWNLOAD/EXTRACT driver ONLY.
#
# This runs the .py steps that fetch and organise GEO Supplementary files:
#   1. validate manifest
#   2. download
#   3. extract archives
#   4. list downloaded files
#
# Everything after this (AnnData loading, inspection, curation, merge, h5ad
# saving) is done in the Jupyter notebooks under notebooks/  -- NOT here.
#
#   ./run.sh                 # steps 00-03
#   ./run.sh 01              # only download
#   PYTHON=../v1/.venv/bin/python ./run.sh
set -euo pipefail
cd "$(dirname "$0")"

if [[ -n "${PYTHON:-}" ]]; then
    PY="$PYTHON"
elif [[ -x ".venv/bin/python" ]]; then
    PY=".venv/bin/python"
else
    PY="python3"
fi

ALL=(00_validate_manifest 01_download_geo_supplement 02_extract_archives \
     03_list_downloaded_files)

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
echo "download/extract done. Continue in notebooks/python/ (start with 00_overview.ipynb)."
