import shutil
from setuptools import setup, Extension
from numpy import get_include as get_numpy_include
from os import path

version = (0,1,0)



def main():
    this_directory = path.abspath(path.dirname(__file__))
    long_description=""
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

    setup(name="rfcnt",
          version="%d.%d.%d" % version,
          description="Python interface for rainflow counting",
          long_description=long_description,
          long_description_content_type='text/markdown',
          keywords='rainflow counting',
          author="Andreas Martin",
          license='BSD-2-Clause License',
          url='http://github.com/AndreasMartin72/rainflow',
          install_requires=['numpy'],
          py_modules=["rfcnt.utils", "rfcnt.examples"],
          packages=["rfcnt", "rfcnt.tests"],
          package_dir={"rfcnt": "", "rfcnt.tests": "tests"},
          package_data={"rfcnt": ["run_examples.py", "run_tests.py", "requirements.txt", "README.md", "LICENSE"],
                        "rfcnt.tests": ["long_series.csv"]},
          ext_modules=[Extension("rfcnt.rfcnt", ["src/rfcnt.cpp", "src/rainflow.c"],
                       define_macros=[
                            ('NPY_NO_DEPRECATED_API',     'NPY_1_7_API_VERSION'),
                            ('RFC_HAVE_CONFIG_H',         '0'),
                            ('RFC_VERSION_MAJOR',         str(version[0])),
                            ('RFC_VERSION_MINOR',         str(version[1])),
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
                            ('RFC_DEBUG_FLAGS',           '0'),
                            ('RFC_EXPORT_MEX',            '0')],
                        include_dirs=['src', get_numpy_include()],
                       )],
          classifiers=[
              'Development Status :: 4 - Beta',
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
          ])

if __name__ == "__main__":
    main()