#!/usr/bin/env python3
"""Verify that every generated, version-carrying file matches CMakeLists.txt.

CMakeLists.txt's `project(rainflow VERSION X.Y.Z ...)` is the single source of
truth for the project version. `cmake -S. -Bbuild` regenerates config.h,
version.py, README.rst and the docs/*.rst files from their .in templates via
configure_file(). This script re-parses CMakeLists.txt the same way and checks
that the checked-in generated files actually match it -- i.e. that someone
bumped the version and forgot to re-run cmake before committing.

Usage:
    python3 tools/check_version_sync.py   # exit 1 and list mismatches, else 0
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CMAKELISTS = ROOT / "CMakeLists.txt"


def parse_expected_version() -> tuple[str, str, str, str, str]:
    text = CMAKELISTS.read_text(encoding="utf-8")

    m = re.search(r"project\s*\(\s*\w+\s+VERSION\s+(\d+)\.(\d+)\.(\d+)", text)
    if not m:
        sys.exit("check_version_sync: could not parse project(VERSION ...) from CMakeLists.txt")
    major, minor, patch = m.group(1), m.group(2), m.group(3)

    m = re.search(r'set\s*\(\s*RFC_VERSION_POSTFIX\s+"([^"]*)"', text)
    postfix = m.group(1) if m else ""

    full = f"{major}.{minor}.{patch}{postfix}"
    return major, minor, patch, postfix, full


def check_config_h(path: Path, major: str, minor: str, patch: str, postfix: str) -> list[str]:
    errors = []
    text = path.read_text(encoding="utf-8")
    checks = [
        ("RFC_VERSION_MAJOR", major),
        ("RFC_VERSION_MINOR", minor),
        ("RFC_VERSION_PATCH", patch),
        ("RFC_VERSION_POSTFIX", postfix),
    ]
    for define, expected in checks:
        m = re.search(rf'#define\s+{define}\s+"([^"]*)"', text)
        if not m:
            errors.append(f"{path}: could not find #define {define}")
        elif m.group(1) != expected:
            errors.append(f"{path}: {define} is \"{m.group(1)}\", expected \"{expected}\"")
    return errors


def check_literal(path: Path, pattern: str, expected: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    matches = re.findall(pattern, text)
    if not matches:
        return [f"{path}: pattern not found: {pattern!r}"]
    errors = []
    for actual in matches:
        if actual != expected:
            errors.append(f"{path}: found version \"{actual}\", expected \"{expected}\"")
    return errors


def main() -> int:
    major, minor, patch, postfix, full = parse_expected_version()

    errors: list[str] = []
    errors += check_config_h(ROOT / "src/lib/config.h", major, minor, patch, postfix)
    errors += check_config_h(ROOT / "src/python/lib/config.h", major, minor, patch, postfix)

    # version.py uses a tuple, not a plain literal -- check its parts individually.
    vpy = ROOT / "src/python/version.py"
    m = re.search(r'_version = \((\d+), (\d+), (\d+), "([^"]*)"\)', vpy.read_text(encoding="utf-8"))
    if not m:
        errors.append(f"{vpy}: could not parse _version tuple")
    elif (m.group(1), m.group(2), m.group(3), m.group(4)) != (major, minor, patch, postfix):
        errors.append(f"{vpy}: _version is {m.groups()}, expected ({major}, {minor}, {patch}, {postfix!r})")

    errors += check_literal(ROOT / "README.rst", r"\.\. \|version\| replace:: (\S+)", full)
    errors += check_literal(ROOT / "docs/index.rst", r"\*\*Documentation Version:\*\* (\S+)", full)
    errors += check_literal(ROOT / "docs/references.rst", r"version = \{([^}]+)\}", full)
    errors += check_literal(ROOT / "docs/installation.rst", r"rfcnt-([^/]+)/rfcnt-\1\.tar\.gz", full)

    if errors:
        print("Version mismatch between CMakeLists.txt and generated files:\n", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        print(
            f"\nCMakeLists.txt declares version {full}. "
            "Run `cmake -S. -Bbuild` to regenerate the files above, then commit the result.",
            file=sys.stderr,
        )
        return 1

    print(f"OK: all generated files match CMakeLists.txt version {full}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
