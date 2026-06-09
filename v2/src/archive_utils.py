"""Safe tar extraction (path-traversal + symlink hardened), recursive nested
tar extraction, and a small file-finder used by the scripts and notebooks.
"""
from __future__ import annotations

import logging
import tarfile
from pathlib import Path

log = logging.getLogger("archive_utils")


# --------------------------------------------------------------------------
# file finding
# --------------------------------------------------------------------------
def find_files(root, patterns=("*",), recursive: bool = True) -> list:
    """Return sorted files under `root` matching any glob in `patterns`."""
    root = Path(root)
    if not root.exists():
        return []
    out: set = set()
    for pat in patterns:
        globber = root.rglob(pat) if recursive else root.glob(pat)
        out.update(p for p in globber if p.is_file())
    return sorted(out)


# --------------------------------------------------------------------------
# safe extraction
# --------------------------------------------------------------------------
def _is_within(directory: Path, target: Path) -> bool:
    directory = directory.resolve()
    target = target.resolve()
    try:
        target.relative_to(directory)
        return True
    except ValueError:
        return False


def _safe_members(tar: tarfile.TarFile, dest: Path) -> list:
    safe = []
    for member in tar.getmembers():
        target = dest / member.name
        if not _is_within(dest, target):
            raise RuntimeError(f"unsafe path in archive: {member.name!r}")
        if member.issym() or member.islnk():
            log.warning("skipping link member %s", member.name)
            continue
        if member.isdev():
            log.warning("skipping device member %s", member.name)
            continue
        safe.append(member)
    return safe


def extract_tar_safe(tar_path, dest) -> Path:
    """Extract one tar safely into `dest`."""
    tar_path, dest = Path(tar_path), Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        members = _safe_members(tar, dest)
        tar.extractall(dest, members=members)
    log.info("extracted %s -> %s (%d members)", tar_path.name, dest, len(members))
    return dest


def is_tar(path) -> bool:
    path = Path(path)
    name = path.name.lower()
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return True
    try:
        return tarfile.is_tarfile(path)
    except Exception:
        return False


def extract_tar_recursive(tar_path, dest, *, max_depth: int = 4) -> Path:
    """Extract a tar then recursively extract any tars found inside it."""
    extract_tar_safe(tar_path, dest)
    if max_depth <= 0:
        return dest
    for inner in list(Path(dest).rglob("*")):
        if inner.is_file() and inner != Path(tar_path) and is_tar(inner):
            inner_dest = inner.parent / (inner.name + "_extracted")
            if inner_dest.exists():
                continue
            log.info("nested archive: %s", inner.name)
            extract_tar_recursive(inner, inner_dest, max_depth=max_depth - 1)
    return dest
