"""GEO supplementary ファイルのダウンロード（レジューム対応・完了済みはスキップ）。

アーカイブ展開は archive_utils.py 側。requests は使わず標準ライブラリ + tqdm のみ。
"""
from __future__ import annotations

import logging
import urllib.error
import urllib.request
from pathlib import Path

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    tqdm = None

log = logging.getLogger("geo_download")

_CHUNK = 1 << 20  # 1 MiB
_TIMEOUT = 120
_UA = "Mozilla/5.0 (compatible; geo-pipeline/1.0)"


def _remote_size(url: str) -> int | None:
    # HEAD で Content-Length を取得（取れなければ None）
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": _UA})
        with urllib.request.urlopen(req, timeout=_TIMEOUT) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:  # pragma: no cover
        return None


def download_file(url: str, dest, *, resume: bool = True, force: bool = False) -> Path:
    """url を dest に保存。`.part` があれば Range で続きから、既存かつサイズ一致ならスキップ。"""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not force:
        remote = _remote_size(url)
        if remote is None or dest.stat().st_size == remote:
            log.info("既存のためスキップ %s", dest.name)
            return dest
        log.warning("サイズ不一致 %s (local=%d remote=%d); 取り直します",
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
        if exc.code == 416:  # Range 不可 = 既に完了
            part.replace(dest)
            return dest
        log.warning("Range 取得に失敗 (%s); %s を最初から取り直し", exc, dest.name)
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
    log.info("ダウンロード完了 %s (%d bytes)", dest.name, dest.stat().st_size)
    return dest


def download_files(file_entries, dest_dir, *, force: bool = False) -> list:
    """manifest の files リストを dest_dir にまとめてダウンロード。"""
    dest_dir = Path(dest_dir)
    out = []
    for entry in file_entries:
        name, url = entry["name"], entry["url"]
        try:
            out.append(download_file(url, dest_dir / name, force=force))
        except Exception as exc:
            if entry.get("optional"):
                log.warning("optional ファイルの失敗を無視して継続: %s (%s)", name, exc)
                continue
            raise RuntimeError(f"ダウンロード失敗 {name} <- {url}: {exc}") from exc
    return out
