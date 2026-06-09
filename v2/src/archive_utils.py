"""tar の安全な展開（パストラバーサル/シンボリックリンク対策）、ネストした
tar の再帰展開、スクリプト/ノートブック共通のファイル探索。
"""
from __future__ import annotations

import logging
import tarfile
from pathlib import Path

log = logging.getLogger("archive_utils")


# --------------------------------------------------------------------------
# ファイル探索
# --------------------------------------------------------------------------
def find_files(root, patterns=("*",), recursive: bool = True) -> list:
    """root 配下で patterns のいずれかに一致するファイルをソートして返す。"""
    root = Path(root)
    if not root.exists():
        return []
    out: set = set()
    for pat in patterns:
        globber = root.rglob(pat) if recursive else root.glob(pat)
        out.update(p for p in globber if p.is_file())
    return sorted(out)


# --------------------------------------------------------------------------
# 安全な展開
# --------------------------------------------------------------------------
def _is_within(directory: Path, target: Path) -> bool:
    # 展開先が dest の外に出ないか確認（パストラバーサル対策）
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
            raise RuntimeError(f"アーカイブ内に不正なパス: {member.name!r}")
        if member.issym() or member.islnk():
            log.warning("リンクメンバをスキップ %s", member.name)
            continue
        if member.isdev():
            log.warning("デバイスメンバをスキップ %s", member.name)
            continue
        safe.append(member)
    return safe


def extract_tar_safe(tar_path, dest) -> Path:
    """tar を1つ安全に dest へ展開する。"""
    tar_path, dest = Path(tar_path), Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        members = _safe_members(tar, dest)
        tar.extractall(dest, members=members)
    log.info("展開 %s -> %s (%d members)", tar_path.name, dest, len(members))
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
    """tar を展開し、中に含まれる tar も再帰的に展開する（GSE178693 用）。"""
    extract_tar_safe(tar_path, dest)
    if max_depth <= 0:
        return dest
    for inner in list(Path(dest).rglob("*")):
        if inner.is_file() and inner != Path(tar_path) and is_tar(inner):
            inner_dest = inner.parent / (inner.name + "_extracted")
            if inner_dest.exists():
                continue
            log.info("ネストしたアーカイブ: %s", inner.name)
            extract_tar_recursive(inner, inner_dest, max_depth=max_depth - 1)
    return dest
