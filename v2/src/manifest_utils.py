"""Manifest loading / validation, project paths, and logging.

Used by both the download scripts (scripts/*.py) and the notebooks.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import yaml

KNOWN_LOADERS = {
    "10x_h5_per_sample", "10x_mtx_per_sample", "combined_umi_tsv_with_metadata",
    "dense_or_text_matrix_bundle", "mtx_or_text_bundle", "dense_gene_by_cell_matrix",
    "processed_count_matrix_with_metadata", "nested_tar_dropseq",
    "R_jupyter_kernel_manual",
}
REQUIRED_KEYS = ["dataset_id", "source_accession", "loader_hint", "output", "files"]


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
def project_root() -> Path:
    """v2/ project root (this file lives in v2/src/)."""
    return Path(__file__).resolve().parent.parent


def project_paths(root: Path | None = None) -> dict:
    root = Path(root) if root is not None else project_root()
    data = root / "data"
    return {
        "root": root,
        "config": root / "config",
        "data": data,
        "raw": data / "raw",
        "extracted": data / "extracted",
        "intermediate_from_r": data / "intermediate_from_r",
        "interim": data / "interim_h5ad",
        "curated": data / "curated_h5ad",
        "merged": data / "merged_h5ad",
        "reports": data / "reports",
    }


def ensure_dirs(paths: dict) -> None:
    for key, value in paths.items():
        if key in ("root", "config"):
            continue
        Path(value).mkdir(parents=True, exist_ok=True)


def dataset_raw_dir(paths: dict, ds: dict) -> Path:
    return Path(paths["raw"]) / ds["source_accession"]


def dataset_extracted_dir(paths: dict, ds: dict) -> Path:
    return Path(paths["extracted"]) / ds["source_accession"]


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------
def load_manifest(path: Path | None = None) -> dict:
    path = Path(path) if path else project_root() / "config" / "dataset_manifest.yaml"
    with open(path) as fh:
        return yaml.safe_load(fh)


def list_datasets(manifest: dict) -> list:
    return manifest.get("datasets", [])


def get_dataset(manifest: dict, key: str) -> dict:
    """Find a dataset by dataset_id, source_accession, or output stem."""
    for ds in manifest.get("datasets", []):
        if key in (ds.get("dataset_id"), ds.get("source_accession")):
            return ds
        if ds.get("output", "").startswith(str(key)):
            return ds
    raise KeyError(f"dataset not found in manifest: {key!r}")


def dataset_files(ds: dict) -> list:
    return ds.get("files", [])


def validate_manifest(manifest: dict) -> list:
    errors: list = []
    seen_ids, seen_out = set(), set()
    for i, ds in enumerate(manifest.get("datasets", [])):
        tag = ds.get("dataset_id", f"<index {i}>")
        for key in REQUIRED_KEYS:
            if not ds.get(key):
                errors.append(f"{tag}: missing required key '{key}'")
        loader = ds.get("loader_hint")
        if loader and loader not in KNOWN_LOADERS:
            errors.append(f"{tag}: unknown loader_hint '{loader}'")
        if ds.get("dataset_id") in seen_ids:
            errors.append(f"{tag}: duplicate dataset_id")
        seen_ids.add(ds.get("dataset_id"))
        if ds.get("output") in seen_out:
            errors.append(f"{tag}: duplicate output")
        seen_out.add(ds.get("output"))
        for f in ds.get("files", []):
            if not f.get("name") or not f.get("url"):
                errors.append(f"{tag}: file entry missing name/url: {f}")
    return errors
