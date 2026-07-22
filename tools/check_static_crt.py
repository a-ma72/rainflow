"""Verify the compiled rfcnt extension does NOT dynamically depend on the
MSVC redistributable runtime DLLs (vcruntime140.dll, vcruntime140_1.dll,
msvcp140.dll).

If it does, the /MT static-CRT compile flag in setup.py's build_ext override
has been dropped or overridden by a later flag -- which reintroduces the
order-dependent kernel crash when this extension shares a process with
another extension built against a different CRT version (e.g. pandas).

Windows-only. No-op (exit 0) on other platforms so it's safe to call
unconditionally from a step that's itself OS-gated in the workflow.
"""
from __future__ import annotations

import sys

FORBIDDEN_DLLS = {"vcruntime140.dll", "vcruntime140_1.dll", "msvcp140.dll"}


def main() -> int:
    if sys.platform != "win32":
        print("Not Windows -- skipping static-CRT check.")
        return 0

    try:
        import pefile
    except ImportError:
        print("pefile is required for this check: pip install pefile", file=sys.stderr)
        return 1

    import rfcnt.rfcnt as ext  # the compiled extension submodule

    pyd_path = ext.__file__
    print(f"Inspecting: {pyd_path}")

    pe = pefile.PE(pyd_path)
    pe.parse_data_directories()

    found = set()
    for entry in getattr(pe, "DIRECTORY_ENTRY_IMPORT", []):
        dll_name = entry.dll.decode("ascii").lower()
        if dll_name in FORBIDDEN_DLLS:
            found.add(dll_name)

    if found:
        print(
            f"FAIL: {pyd_path} dynamically imports {sorted(found)}.\n"
            f"The /MT static-CRT flag appears to have been dropped from "
            f"setup.py's msvc build_ext flags, or a later flag is overriding it.\n"
            f"This will reintroduce the pandas/rfcnt import-order kernel crash.",
            file=sys.stderr,
        )
        return 1

    print("OK: no dynamic dependency on vcruntime140/msvcp140 -- static CRT confirmed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
