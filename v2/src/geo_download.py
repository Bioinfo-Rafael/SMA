"""GEO supplementary download + safe archive extraction.

* Downloads with HTTP range-resume; already-complete files are skipped.
* tar extraction is hardened against path traversal and symlink escapes.
* Nested tars (e.g. GSE178693 RAW.tar containing per-sample tars) are
  extracted recursively.

Only the Python standard library + tqdm are used (no `requests` dependency).
"""
from __future__ import annotations

import logging
import os
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover - tqdm optional at import time
    tqdm = None

log = logging.getLogger("geo_download")

_CHUNK = 1 << 20  # 1 MiB
_TIMEOUT = 120
_UA = "Mozilla/5.0 (compatible; geo-pipeline/1.0)"


# --------------------------------------------------------------------------
# download
# --------------------------------------------------------------------------
def _remote_size(url: str) -> int | None:
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": _UA})
        with urllib.request.urlopen(req, timeout=_TIMEOUT) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:  # pragma: no cover - HEAD not always allowed
        return None


def download_file(url: str, dest, *, resume: bool = True, force: bool = False) -> Path:
    """Download `url` to `dest`. Resumes a partial `.part` file when possible,
    skips when `dest` already exists (unless force)."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not force:
        remote = _remote_size(url)
        if remote is None or dest.stat().st_size == remote:
            log.info("skip existing %s", dest.name)
            return dest
        log.warning("size mismatch for %s (local=%d remote=%d); re-downloading",
                    dest.name, dest.stat().st_size, remote)

    part = dest.with_name(dest.name + ".part")
    existing = part.stat().st_size if (resume and part.exists()) else 0

    headers = {"User-Agent": _UA}
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"

    try:
        resp = urllib.request.urlopen(
            urllib.request.Request(url, headers=headers), timeout=_TIMEOUT)
    except urllib.error.HTTPError as exc:
        if exc.code == 416:  # range not satisfiable -> already complete
            part.replace(dest)
            return dest
        log.warning("range request failed (%s); restarting %s", exc, dest.name)
        existing = 0
        if part.exists():
            part.unlink()
        resp = urllib.request.urlopen(
            urllib.request.Request(url, headers={"User-Agent": _UA}), timeout=_TIMEOUT)

    resumed = existing > 0 and getattr(resp, "status", 200) == 206
    mode = "ab" if resumed else "wb"
    if not resumed:
        existing = 0

    total = resp.headers.get("Content-Length")
    total = (int(total) + existing) if total is not None else None

    bar = None
    if tqdm is not None:
        bar = tqdm(total=total, initial=existing, unit="B", unit_scale=True,
                   desc=dest.name, leave=False)
    try:
        with open(part, mode) as fh:
            while True:
                chunk = resp.read(_CHUNK)
                if not chunk:
                    break
                fh.write(chunk)
                if bar is not None:
                    bar.update(len(chunk))
    finally:
        if bar is not None:
            bar.close()
        resp.close()

    part.replace(dest)
    log.info("downloaded %s (%d bytes)", dest.name, dest.stat().st_size)
    return dest


def download_files(file_entries, dest_dir, *, force: bool = False) -> list[Path]:
    """Download a list of manifest file entries into dest_dir."""
    dest_dir = Path(dest_dir)
    out: list[Path] = []
    for entry in file_entries:
        name = entry["name"]
        url = entry["url"]
        try:
            out.append(download_file(url, dest_dir / name, force=force))
        except Exception as exc:
            if entry.get("optional"):
                log.warning("optional file failed, continuing: %s (%s)", name, exc)
                continue
            raise RuntimeError(f"download failed for {name} <- {url}: {exc}") from exc
    return out


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


def _safe_members(tar: tarfile.TarFile, dest: Path):
    safe = []
    for member in tar.getmembers():
        target = dest / member.name
        if not _is_within(dest, target):
            raise RuntimeError(f"unsafe path in archive: {member.name!r}")
        if member.issym() or member.islnk():
            log.warning("skipping link member %s in %s", member.name, tar.name)
            continue
        if member.isdev():
            log.warning("skipping device member %s", member.name)
            continue
        safe.append(member)
    return safe


def safe_extract_tar(tar_path, dest) -> Path:
    """Extract one tar safely (path-traversal + link hardened)."""
    tar_path = Path(tar_path)
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        members = _safe_members(tar, dest)
        tar.extractall(dest, members=members)
    log.info("extracted %s -> %s (%d members)", tar_path.name, dest, len(members))
    return dest


def _is_tar(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return True
    try:
        return tarfile.is_tarfile(path)
    except Exception:
        return False


def safe_extract_recursive(tar_path, dest, *, max_depth: int = 4) -> Path:
    """Extract a tar, then recursively extract any tars found inside it
    (each into a sibling directory named after the inner archive)."""
    safe_extract_tar(tar_path, dest)
    if max_depth <= 0:
        return dest
    for inner in list(Path(dest).rglob("*")):
        if inner.is_file() and inner != Path(tar_path) and _is_tar(inner):
            inner_dest = inner.parent / (inner.name + "_extracted")
            if inner_dest.exists():
                continue
            log.info("nested archive: %s", inner.name)
            safe_extract_recursive(inner, inner_dest, max_depth=max_depth - 1)
    return dest
