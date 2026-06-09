#!/usr/bin/env bash
# v2 ダウンロード/展開 専用ドライバ。
#
# GEO Supplementary を取得・整理する .py ステップだけを実行する:
#   1. manifest 検証
#   2. ダウンロード
#   3. アーカイブ展開
#   4. ダウンロード済みファイル一覧
#   5. 俯瞰（overview）
#
# これ以降（AnnData 化・確認・curate・merge・h5ad 保存）は notebooks/ の
# Jupyter ノートブックで行う。ここでは実行しない。
#
#   ./run.sh                 # ステップ 00-04
#   ./run.sh 01              # ダウンロードだけ
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
     03_list_downloaded_files 04_overview)

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
echo "ダウンロード/展開 完了。続きは notebooks/python/01_load_and_inspect_each_gse.ipynb から。"
