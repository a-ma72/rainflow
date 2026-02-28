# Source distribution (./dist)
# python setup.py build sdist
#
# Binary distribution (./build)
# python -m build -nwx
# python setup.py bdist_wheel --plat-name=win-amd64
# python setup.py bdist --formats=wininst
# python setup.py bdist_wininst --title= --bitmap=
# pip install --force-reinstall --no-deps --no-build-isolation package.tar
# python setup.py build_clib -f
# python setup.py build_ext -fi  # Build pyd modules for folder `_ext`
# python -mbuild -n

# If build with mingw32 compiler (TDM-GCC64):
# Comment out the get_msvcr() occurrences in PATHONPATH/Lib/distutils/cygwinccompiler.py
# Create a file distutils.cfg in PYTHON_ROOT/Lib/distutils:
# [build]
# compiler=mingw32

# PyPi
# python3 -m pip install --upgrade build twine
# python3 -m build
# python3 -m twine upload --repository testpypi dist/*
# pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple rfcnt==0.2.0
# python3 -m twine upload --repository pypi dist/rfcnt-0.5.1rc1*


# Install (Jupyter Notebook)
# !export CFLAGS='-std=c++11' && pip install rfcnt

from __future__ import annotations

import logging
import re
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext

try:
    # Check if numpy is installed
    from numpy import __version__ as np_version
    from numpy import get_include as np_get_include
except ImportError:
    np_version = "NUMPY_NOTFOUND"
    def np_get_include() -> str:
        """Return a placeholder string indicating that NumPy is not installed."""
        return "NUMPY_NOTFOUND"

class build_ext(_build_ext):  # noqa: N801
    """Custom build_ext that applies per-file/per-language compiler flags."""

    def build_extension(self, ext: Extension) -> None:
        ct = self.compiler.compiler_type
        
        # Original logic: modify the extension object before building
        for i, src in enumerate(ext.sources):
            # Check file extension to determine flags
            is_cpp = src.endswith((".cpp", ".cc", ".cxx"))
            
            if ct == "msvc":
                if is_cpp:
                    ext.extra_compile_args.append("/std:c++14")
                ext.extra_compile_args.extend([
                    "/wd4100", "/wd4101", "/wd4189", "/wd4505"
                ])
            elif ct in ("unix", "mingw32"):
                # Use -Wno-... to SUPPRESS warnings
                ext.extra_compile_args.extend([
                    "-Wno-unused-variable",
                    "-Wno-unused-function",
                    "-Wno-unused-but-set-variable",
                    "-Wno-switch",
                    "-Wno-unused-value",
                ])
                
                # IMPORTANT: Only add -std=c99 to .c files and -std=c++11 to .cpp files
                # Note: On many Clang versions, we skip adding -std=c99 to the global 
                # extra_compile_args if C++ files are present. 
                # Better: Use the compiler's ability to handle this by file type.
                if is_cpp and "-std=c++11" not in ext.extra_compile_args:
                    ext.extra_compile_args.append("-std=c++11")
                elif not is_cpp and "-std=c99" not in ext.extra_compile_args:
                    # To avoid the Clang error, we only add this if the extension 
                    # is purely C, OR we rely on the fact that modern Clang 
                    # defaults to a compatible C standard. 
                    # If you MUST have c99 for the .c file:
                    ext.extra_compile_args.append("-std=c99")

        # If Clang is still complaining because of the mix, we remove -std=c99 
        # and let the compiler use its default, which is usually newer/compatible.
        if ct in ("unix", "mingw32") and any(s.endswith(".cpp") for s in ext.sources):
            if "-std=c99" in ext.extra_compile_args:
                ext.extra_compile_args.remove("-std=c99")

        super().build_extension(ext)


def parse_version_file(version_file_path: Path) -> tuple[str, str, str]:
    """Parse version.py and return _version, __version__, and __author__."""
    namespace: dict = {}
    exec(version_file_path.read_text(), namespace)  # noqa: S102
    _version = namespace.get("_version")
    __version__ = namespace.get("__version__")
    __author__ = namespace.get("__author__")
    if _version is None or __version__ is None or __author__ is None:
        msg = "Could not parse version or author information from version.py"
        raise ValueError(msg)
    return _version, __version__, __author__


def main() -> None:
    """Set up the rfcnt package by reading metadata, version info, and configuring the build process."""
    this_directory = Path(__file__).parent.resolve()
    long_description = ""
    with (this_directory / "README.rst").open(encoding="utf-8") as f:
        long_description = f.read()
    with (this_directory / "lib" / "config.h").open() as f:
        RFC_VERSION_MAJOR = RFC_VERSION_MINOR = None
        for line in f:
            match = re.match(r".*#define\s+RFC_VERSION_(MAJOR|MINOR)\s+\"(\d+)\".*", line)
            if match:
                if match.group(1) == "MAJOR":
                    RFC_VERSION_MAJOR = match.group(2)
                else:
                    RFC_VERSION_MINOR = match.group(2)
        if RFC_VERSION_MAJOR is None:
            msg = "Can't locate version signature: RFC_VERSION_MAJOR is missing."
            raise ValueError(msg)
        if RFC_VERSION_MINOR is None:
            msg = "Can't locate version signature: RFC_VERSION_MINOR is missing."
            raise ValueError(msg)
    version_file = this_directory / "version.py"
    _version, __version__, __author__ = parse_version_file(version_file)

    define_macros = [
        ("NPY_NO_DEPRECATED_API",     "NPY_1_7_API_VERSION"),
        ("NPY_TARGET_VERSION",        "NPY_1_19_API_VERSION"),
        ("RFC_HAVE_CONFIG_H",         "0"),
        ("RFC_VERSION_MAJOR",         RFC_VERSION_MAJOR),
        ("RFC_VERSION_MINOR",         RFC_VERSION_MINOR),
        ("RFC_USE_INTEGRAL_COUNTS",   "0"),
        ("RFC_USE_HYSTERESIS_FILTER", "1"),
        ("RFC_MINIMAL",               "0"),
        ("RFC_TP_SUPPORT",            "1"),
        ("RFC_HCM_SUPPORT",           "1"),
        ("RFC_ASTM_SUPPORT",          "1"),
        ("RFC_USE_DELEGATES",         "1"),
        ("RFC_GLOBAL_EXTREMA",        "1"),
        ("RFC_DAMAGE_FAST",           "1"),
        ("RFC_DH_SUPPORT",            "1"),
        ("RFC_AT_SUPPORT",            "1"),
        ("RFC_AR_SUPPORT",            "1"),
        ("RFC_DEBUG_FLAGS",           "0"),
        ("RFC_EXPORT_MEX",            "0"),
        ("RFC_EXPORT_PY",             "1")]

    setup(
        name="rfcnt",
        version=__version__,
        description="Python interface for rainflow counting",
        long_description=long_description,
        long_description_content_type="text/markdown",
        keywords="rainflow counting",
        author=__author__,
        license="BSD-2-Clause",
        url="http://github.com/a-ma72/rainflow",
        install_requires=["numpy"],
        packages=["rfcnt", "rfcnt.tests"],
        package_dir={"rfcnt": "", "rfcnt.tests": "tests"},
        package_data={
            "rfcnt": ["*.py*", "requirements.txt", "README.rst", "LICENSE"],
            "rfcnt.tests": ["*.py", "long_series.csv"],
        },
        cmdclass={"build_ext": build_ext},
        ext_modules=[
            Extension(
                "rfcnt.rfcnt", ["src/rfcnt.cpp", "lib/rainflow.c"],
                define_macros=define_macros,
                include_dirs=["src", "lib", np_get_include()],
            ),
        ],
        classifiers=[
            "Development Status :: 6 - Mature",
            "Environment :: Console",
            "Framework :: Buildout :: Extension",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "Intended Audience :: Information Technology",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX",
            "Programming Language :: Python :: 3",
            "Programming Language :: C++",
            "Programming Language :: C",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Physics",
        ],
    )



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    main()
    logger.info("numpy: %s", np_version)
