"""Display helpers for the inspection notebook (notebooks/python/02_*).

These print to the notebook and return DataFrames so you can sort/filter. They
deliberately do NOT pick cell-type / cluster / condition columns -- that is the
human's job after looking at the output.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

DEFAULT_SAMPLE_KEYS = [
    "sample_id", "source_accession", "disease_status", "treatment",
    "tissue", "enrichment", "data_status",
]


def summarize_adata(adata, name: str | None = None) -> dict:
    """Print a compact summary and return it as a dict."""
    X = adata.X
    issparse = sp.issparse(X)
    info = {
        "name": name,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "X_sparse": issparse,
        "X_dtype": str(getattr(X, "dtype", "NA")),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "layers": list(adata.layers.keys()),
        "obsm": list(adata.obsm.keys()),
        "uns_keys": list(adata.uns.keys()),
    }
    try:
        info["X_min"] = float(X.min()) if issparse else float(np.min(X))
        info["X_max"] = float(X.max()) if issparse else float(np.max(X))
        info["X_sum"] = float(X.sum())
    except Exception as exc:  # pragma: no cover
        info["X_stats_error"] = str(exc)

    head = f"=== {name} ===\n" if name else ""
    print(f"{head}{adata.n_obs} cells x {adata.n_vars} genes | "
          f"X sparse={info['X_sparse']} dtype={info['X_dtype']} "
          f"min={info.get('X_min')} max={info.get('X_max')} sum={info.get('X_sum')}")
    print(f"obs columns ({len(info['obs_columns'])}): {info['obs_columns']}")
    print(f"var columns ({len(info['var_columns'])}): {info['var_columns']}")
    print(f"layers={info['layers']} obsm={info['obsm']} uns={info['uns_keys']}")
    return info


def show_obs_columns(adata) -> pd.DataFrame:
    rows = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        rows.append({"column": col, "dtype": str(s.dtype),
                     "non_null": int(s.notna().sum()),
                     "n_missing": int(s.isna().sum()),
                     "n_unique": int(s.nunique(dropna=True)),
                     "example": str(s.iloc[0]) if len(s) else ""})
    return pd.DataFrame(rows)


def show_var_columns(adata) -> pd.DataFrame:
    rows = []
    for col in adata.var.columns:
        s = adata.var[col]
        rows.append({"column": col, "dtype": str(s.dtype),
                     "non_null": int(s.notna().sum()),
                     "n_unique": int(s.nunique(dropna=True)),
                     "example": str(s.iloc[0]) if len(s) else ""})
    return pd.DataFrame(rows)


def show_obs_value_counts(adata, max_unique: int = 50, top: int = 30) -> dict:
    """Print value_counts for low-cardinality / categorical / object obs columns.
    Returns {column: value_counts Series}."""
    out: dict = {}
    for col in adata.obs.columns:
        s = adata.obs[col]
        is_cat = str(s.dtype) == "category"
        if s.dtype.kind in ("O", "b") or is_cat or (s.nunique(dropna=True) <= max_unique):
            vc = s.astype(str).value_counts().head(top)
            out[col] = vc
            print(f"\n--- {col} ({s.nunique(dropna=True)} unique) ---")
            print(vc.to_string())
    return out


def show_numeric_obs_summary(adata) -> pd.DataFrame:
    num = adata.obs.select_dtypes(include=[np.number])
    if num.shape[1] == 0:
        print("(no numeric obs columns)")
        return pd.DataFrame()
    desc = num.describe().T
    print(desc.to_string())
    return desc


def show_sample_counts(adata, keys=None) -> dict:
    """Print cell counts per value for each key present in obs."""
    keys = keys or DEFAULT_SAMPLE_KEYS
    out: dict = {}
    for key in keys:
        if key in adata.obs.columns:
            vc = adata.obs[key].astype(str).value_counts()
            out[key] = vc
            print(f"\n=== cells by {key} ===")
            print(vc.to_string())
    return out
