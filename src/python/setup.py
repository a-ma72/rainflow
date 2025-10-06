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
# python3 -m twine upload --repository pypi dist/rfcnt-0.5.0rc1*


# Install (Jupyter Notebook)
# !export CFLAGS='-std=c++11' && pip install rfcnt

import re
from os import path
from setuptools import setup, Extension


try:
    # Check if numpy is installed
    from numpy import __version__ as np_version
    from numpy import get_include as np_get_include
except ImportError:
    np_version = "NUMPY_NOTFOUND"
    def np_get_include():
        return "NUMPY_NOTFOUND"


def main():
    this_directory = path.abspath(path.dirname(__file__))
    long_description = ""
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
    with open(path.join(this_directory, 'lib', 'config.h')) as f:
        RFC_VERSION_MAJOR = RFC_VERSION_MINOR = None
        for line in f:
            match = re.match(r".*#define\s+RFC_VERSION_(MAJOR|MINOR)\s+\"(\d+)\".*", line)
            if match:
                if match.group(1) == "MAJOR":
                    RFC_VERSION_MAJOR = match.group(1)
                else:
                    RFC_VERSION_MINOR = match.group(2)
        assert RFC_VERSION_MAJOR is not None and RFC_VERSION_MINOR is not None, \
            "Can't locate version signature."
    ver_ns = {}
    with open(path.join(this_directory, "version.py")) as f:
        exec(f.read(), ver_ns)
    _version = ver_ns["_version"]
    __version__ = ver_ns["__version__"]
    __author__ = ver_ns["__author__"]
    del ver_ns, f

    define_macros = [
        ('NPY_NO_DEPRECATED_API',     'NPY_1_7_API_VERSION'),
        ('RFC_HAVE_CONFIG_H',         '0'),
        ('RFC_VERSION_MAJOR',         RFC_VERSION_MAJOR),
        ('RFC_VERSION_MINOR',         RFC_VERSION_MINOR),
        ('RFC_USE_INTEGRAL_COUNTS',   '0'),
        ('RFC_USE_HYSTERESIS_FILTER', '1'),
        ('RFC_MINIMAL',               '0'),
        ('RFC_TP_SUPPORT',            '1'),
        ('RFC_HCM_SUPPORT',           '1'),
        ('RFC_ASTM_SUPPORT',          '1'),
        ('RFC_USE_DELEGATES',         '1'),
        ('RFC_GLOBAL_EXTREMA',        '1'),
        ('RFC_DAMAGE_FAST',           '1'),
        ('RFC_DH_SUPPORT',            '1'),
        ('RFC_AT_SUPPORT',            '1'),
        ('RFC_AR_SUPPORT',            '1'),
        ('RFC_DEBUG_FLAGS',           '0'),
        ('RFC_EXPORT_MEX',            '0'),
        ('RFC_EXPORT_PY',             '0')]

    setup(
        name="rfcnt",
        version=__version__,
        description="Python interface for rainflow counting",
        long_description=long_description,
        long_description_content_type='text/markdown',
        keywords='rainflow counting',
        author=__author__,
        license='BSD-2-Clause',
        url='http://github.com/a-ma72/rainflow',
        setup_requires=['wheel'],
        install_requires=['numpy'],
        packages=["rfcnt", "rfcnt.tests"],
        package_dir={"rfcnt": "", "rfcnt.tests": "tests"},
        package_data={
            "rfcnt": ["*.py*", "_ext/*", "prebuilds.json",
                      "requirements.txt", "README.md", "LICENSE"],
            "rfcnt.tests": ["*.py", "long_series.csv"]
        },
        libraries=[("rainflow_c", {"sources": ["lib/rainflow.c"],
                                   "macros": define_macros})],
        ext_modules=[
            Extension(
                "rfcnt.rfcnt", ["src/rfcnt.cpp"],
                define_macros=define_macros,
                include_dirs=['src', 'lib', np_get_include()],
                extra_compile_args=['-std=c++11'],
            )
        ],
        classifiers=[
            'Development Status :: 6 - Mature',
            'Environment :: Console',
            'Framework :: Buildout :: Extension',
            'Intended Audience :: Developers',
            'Intended Audience :: Education',
            'Intended Audience :: Information Technology',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 3',
            'Programming Language :: C++',
            'Programming Language :: C',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Information Analysis',
            'Topic :: Scientific/Engineering :: Physics',
        ]
    )


if __name__ == "__main__":
    main()
    print("numpy: %s" % (np_version,))

