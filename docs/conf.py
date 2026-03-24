"""Sphinx configuration — version is read from the root CMakeLists.txt."""

import re
from pathlib import Path

# Read version from CMakeLists.txt (single source of truth)
_cmake = (Path(__file__).parent.parent / "CMakeLists.txt").read_text()
_match = re.search(r"project\s*\(\s*\w+\s+VERSION\s+([\d.]+)", _cmake)
if not _match:
    raise RuntimeError("Could not parse version from CMakeLists.txt")
_postfix_match = re.search(r'set\s*\(\s*RFC_VERSION_POSTFIX\s+"([^"]*)"', _cmake)
_postfix = _postfix_match.group(1) if _postfix_match else ""
version = release = _match.group(1) + _postfix

project = "rainflow"
author = "Andreas Martin"
copyright = f"2026, {author}"

extensions = []
master_doc = "index"

# Make |version| and |release| available in all RST files without explicit import
rst_prolog = f"""
.. |version| replace:: {version}
.. |release| replace:: {release}
"""
