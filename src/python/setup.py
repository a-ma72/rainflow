import logging
import re
from pathlib import Path

from setuptools import Extension, setup

try:
    # Check if numpy is installed
    from numpy import __version__ as np_version
    from numpy import get_include as np_get_include
except ImportError:
    np_version = "NUMPY_NOTFOUND"
    def np_get_include() -> str:
        """Return a placeholder string indicating that NumPy is not installed."""
        return "NUMPY_NOTFOUND"


def parse_version_file(version_file_path: Path) -> tuple[str, str, str]:
    """Parse version.py and return _version, __version__, and __author__."""
    _version = __version__ = __author__ = None
    with version_file_path.open() as f:
        for line in f:
            if line.startswith("_version"):
                _version = line.split("=")[1].strip().strip("'\"")
            elif line.startswith("__version__"):
                __version__ = line.split("=")[1].strip().strip("'\"")
            elif line.startswith("__author__"):
                __author__ = line.split("=")[1].strip().strip("'\"")
    if _version is None or __version__ is None or __author__ is None:
        msg = "Could not parse version or author information from version.py"
        raise ValueError(msg)
    return _version, __version__, __author__


def main() -> None:
    """Set up the rfcnt package by reading metadata, version info, and configuring the build process."""
    this_directory = Path(__file__).parent.resolve()
    long_description = ""
    with (this_directory / "README.md").open(encoding="utf-8") as f:
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
        ("RFC_EXPORT_PY",             "0")]

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
        setup_requires=["wheel"],
        install_requires=["numpy"],
        packages=["rfcnt", "rfcnt.tests"],
        package_dir={"rfcnt": "", "rfcnt.tests": "tests"},
        package_data={
            "rfcnt": ["*.py*", "requirements.txt", "README.md", "LICENSE"],
            "rfcnt.tests": ["*.py", "long_series.csv"],
        },
        ext_modules=[
            Extension(
                "rfcnt.rfcnt", ["src/rfcnt.cpp", "lib/rainflow.c"],
                define_macros=define_macros,
                include_dirs=["src", "lib", np_get_include()],
                extra_compile_args=["-std=c++11"],
                # extra_compile_args={
                #     'c': ['-O2'],          # for C files
                #     'c++': ['-std=c++11'], # for C++ files
                # },
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
