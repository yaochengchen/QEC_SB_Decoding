from __future__ import annotations

import hashlib
import platform
from pathlib import Path
import subprocess
import sys


def dependency_versions() -> dict:
    import matplotlib
    import networkx
    import numpy
    import pandas
    import pymatching
    import scipy
    import stim

    return {
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "matplotlib": matplotlib.__version__,
        "networkx": networkx.__version__,
        "stim": stim.__version__,
        "pymatching": pymatching.__version__,
        "platform": platform.platform(),
    }


def source_manifest(root: Path | None = None) -> dict:
    root = Path(root or Path(__file__).resolve().parents[2]).resolve()
    files = sorted(
        path for folder in (root / "qec_sb", root / "scripts")
        for path in folder.rglob("*.py")
        if "__pycache__" not in path.parts
    )
    entries = []
    combined = hashlib.sha256()
    for path in files:
        relative = path.relative_to(root).as_posix()
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        entries.append({"path": relative, "sha256": digest})
        combined.update(relative.encode("utf-8"))
        combined.update(b"\0")
        combined.update(digest.encode("ascii"))
        combined.update(b"\n")
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    return {
        "git_revision": revision,
        "git_dirty": dirty,
        "source_manifest_sha256": combined.hexdigest(),
        "files": entries,
    }
