"""Logging setup and the manual-inspection report writer.

The inspection report exposes *every* obs/var column, dtype, missing count,
value counts and numeric summaries so a human can decide cell-type / cluster
columns later. It never picks those columns itself.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def get_logger(name: str = "pipeline", level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s", "%H:%M:%S"))
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.propagate = False
    return logger


# columns whose per-group cell counts are always reported when present
GROUP_KEYS = [
    "sample_id", "source_accession", "disease_status", "treatment", "enrichment",
]


def inspect_anndata(adata, dataset_id: str, reports_dir) -> dict:
    """Write the 4 report artefacts for one interim h5ad and return their paths."""
    reports_dir = Path(reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)
    vc_dir = reports_dir / f"{dataset_id}_obs_value_counts"
    vc_dir.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []

    def emit(text: str = "") -> None:
        lines.append(str(text))

    emit(f"dataset_id : {dataset_id}")
    emit(f"n_obs      : {adata.n_obs}")
    emit(f"n_vars     : {adata.n_vars}")
    emit("")
    emit(f"obs_columns: {list(adata.obs.columns)}")
    emit(f"var_columns: {list(adata.var.columns)}")
    emit(f"layers     : {list(adata.layers.keys())}")
    emit(f"obsm       : {list(adata.obsm.keys())}")
    emit(f"uns_keys   : {list(adata.uns.keys())}")
    emit("")

    X = adata.X
    issparse = sp.issparse(X)
    emit(f"X_sparse   : {issparse}")
    emit(f"X_dtype    : {getattr(X, 'dtype', 'NA')}")
    try:
        xmin = X.min() if issparse else np.min(X)
        xmax = X.max() if issparse else np.max(X)
        xsum = X.sum()
        emit(f"X_min      : {xmin}")
        emit(f"X_max      : {xmax}")
        emit(f"X_sum      : {xsum}")
    except Exception as exc:  # pragma: no cover
        emit(f"X_stats    : ERROR {exc}")
    emit("")

    # ---- obs columns: dtype / non-null / uniques / value counts / describe ----
    obs_rows = []
    emit("=== obs columns ===")
    for col in adata.obs.columns:
        series = adata.obs[col]
        non_null = int(series.notna().sum())
        n_unique = int(series.nunique(dropna=True))
        dtype = str(series.dtype)
        obs_rows.append({"column": col, "dtype": dtype,
                         "non_null": non_null, "n_missing": int(series.isna().sum()),
                         "n_unique": n_unique})
        emit(f"  {col:24s} dtype={dtype:12s} non_null={non_null:>8d} "
             f"unique={n_unique:>6d}")

        kind = series.dtype.kind
        is_categorical = str(series.dtype) == "category"
        if kind in ("O", "b") or is_categorical or kind in ("i", "u"):
            vc = series.astype(str).value_counts().head(30)
            vc.rename("count").to_csv(vc_dir / f"{col}.csv", header=True)
        if kind in ("f", "i", "u"):
            try:
                desc = series.describe().to_dict()
                emit(f"      describe: {desc}")
            except Exception:  # pragma: no cover
                pass
    pd.DataFrame(obs_rows).to_csv(
        reports_dir / f"{dataset_id}_obs_columns.csv", index=False)
    emit("")

    # ---- var columns ----
    var_rows = []
    emit("=== var columns ===")
    for col in adata.var.columns:
        series = adata.var[col]
        var_rows.append({"column": col, "dtype": str(series.dtype),
                         "non_null": int(series.notna().sum()),
                         "n_unique": int(series.nunique(dropna=True))})
        emit(f"  {col:24s} dtype={str(series.dtype):12s} "
             f"non_null={int(series.notna().sum()):>8d}")
    pd.DataFrame(var_rows).to_csv(
        reports_dir / f"{dataset_id}_var_columns.csv", index=False)
    emit("")

    # ---- grouped cell counts ----
    for key in GROUP_KEYS:
        if key in adata.obs.columns:
            emit(f"=== cells by {key} ===")
            for value, count in adata.obs[key].astype(str).value_counts().items():
                emit(f"  {value}: {count}")
            emit("")

    summary_path = reports_dir / f"{dataset_id}_summary.txt"
    summary_path.write_text("\n".join(lines))
    return {
        "summary": summary_path,
        "obs_columns": reports_dir / f"{dataset_id}_obs_columns.csv",
        "var_columns": reports_dir / f"{dataset_id}_var_columns.csv",
        "value_counts_dir": vc_dir,
    }
