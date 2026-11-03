"""Repository integrity: every nanofab module imported by tracked code is itself tracked by git.

A module that exists only in a developer's working tree imports fine locally but breaks a clean
checkout (src/nanofab/whole_wick/ was imported by committed code before it was committed).
"""

import ast
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from nanofab.provenance import REPO_ROOT

SRC = REPO_ROOT / "src"


def _git_tracked() -> set[str]:
    if shutil.which("git") is None:
        pytest.skip("git not available")
    proc = subprocess.run(["git", "-C", str(REPO_ROOT), "ls-files", "-z"], capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.skip("not a git work tree")
    return {p for p in proc.stdout.split("\0") if p}


def _module_of(path: Path) -> tuple[str, bool] | None:
    """Dotted module name for a file under src/, and whether it is a package __init__."""
    try:
        rel = path.relative_to(SRC)
    except ValueError:
        return None
    parts = list(rel.with_suffix("").parts)
    is_pkg = parts[-1] == "__init__"
    if is_pkg:
        parts.pop()
    return ".".join(parts), is_pkg


def _imported_modules(path: Path) -> set[str]:
    """Absolute dotted names of nanofab modules (and possible submodules) imported by ``path``."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    here = _module_of(path)
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if here is None:
                    continue
                name, is_pkg = here
                pkg = name.split(".") if is_pkg else name.split(".")[:-1]
                base = pkg[: len(pkg) - (node.level - 1)]
                mod = ".".join(base + ([node.module] if node.module else []))
            else:
                mod = node.module or ""
            found.add(mod)
            found.update(f"{mod}.{a.name}" for a in node.names if a.name != "*")
    return {m for m in found if m == "nanofab" or m.startswith("nanofab.")}


def _exists(path: Path) -> bool:
    """Case-sensitive existence check (``from pkg import Grid`` must not match grid.py on Windows)."""
    return path.parent.is_dir() and path.name in os.listdir(path.parent)


def _files_for(module: str) -> list[Path] | None:
    """Source files a module import executes: each parent package __init__ plus the module itself."""
    parts = module.split(".")
    files = []
    for i in range(1, len(parts) + 1):
        stem = SRC.joinpath(*parts[:i])
        if _exists(stem) and _exists(stem / "__init__.py"):
            files.append(stem / "__init__.py")
        elif _exists(stem.with_suffix(".py")):
            files.append(stem.with_suffix(".py"))
        else:
            return None
    return files


def test_imported_nanofab_modules_are_tracked():
    tracked = _git_tracked()
    problems = []
    for rel in sorted(t for t in tracked if t.endswith(".py")):
        path = REPO_ROOT / rel
        if not path.exists():
            continue
        for module in sorted(_imported_modules(path)):
            files = _files_for(module)
            if files is None:
                # ``from pkg import name`` also yields pkg.name, which is often an attribute.
                if _files_for(module.rsplit(".", 1)[0]) is None:
                    problems.append(f"{rel}: imports {module}, which does not exist")
                continue
            for f in files:
                frel = f.relative_to(REPO_ROOT).as_posix()
                if frel not in tracked:
                    problems.append(f"{rel}: imports {module}, but {frel} is not tracked by git")
    assert not problems, "\n".join(sorted(set(problems)))
