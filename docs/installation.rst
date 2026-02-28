============
Installation
============

This page describes how to build and install the rainflow counting package
from sources for different target environments.

Prerequisites
=============

Required Tools
--------------

- **CMake** (3.10 or higher)
- **C/C++ Compiler** with C99 support

  - GCC 4.8+
  - Clang 3.4+
  - MSVC 2015+

Optional Dependencies
---------------------

For Python support:

- **Python** 3.6 or higher
- **NumPy** (via ``oldest-supported-numpy``)
- **setuptools**, **wheel**, **build** packages

For MATLAB support:

- **MATLAB** R2017b or higher

Building from Sources
=====================

Full Build (All Components)
----------------------------

Build all components including C library, Python extension, MATLAB MEX,
and unit tests:

.. code-block:: bash

   cmake -S. -Bbuild -G "Visual Studio 16 2019"
   cmake --build build --config Release

On Linux/macOS:

.. code-block:: bash

   cmake -S. -Bbuild
   cmake --build build --config Release

MATLAB Integration
==================

If you want to use a specific MATLAB installation, set the ``Matlab_ROOT_DIR``
environment variable before running CMake.

Using Environment Variable
---------------------------

.. code-block:: bash

   # Linux/macOS
   export Matlab_ROOT_DIR=/usr/local/MATLAB/R2019b
   cmake -S. -Bbuild
   cmake --build build --config Release

   # Windows PowerShell
   $env:Matlab_ROOT_DIR="C:\Program Files\MATLAB\R2019b"
   cmake -S. -Bbuild -G "Visual Studio 16 2019"
   cmake --build build --config Release

Using CMake Option
------------------

Alternatively, pass it directly to CMake:

.. code-block:: bash

   cmake -S. -Bbuild -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2017b
   cmake --build build --config Release

MATLAB Only Build
-----------------

To build only the MATLAB MEX file:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_UNIT_TEST=0 -G "Visual Studio 16 2019"
   cmake --build build --target rfc_mex --config Release

Python Integration
==================

Python Extension via CMake
---------------------------

To build only the Python extension using CMake:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_EXPORT_MEX=0 -DRFC_UNIT_TEST=0
   cmake --build build --target rfcnt --config Release

Python Wheel Package
--------------------

To build a Python wheel package (recommended for Python users):

.. code-block:: bash

   cd src/python
   pip install setuptools build wheel oldest-supported-numpy
   python -m build -nw

This creates a wheel file in ``src/python/dist/`` that can be installed with pip:

.. code-block:: bash

   pip install dist/rfcnt-*.whl

Installing from PyPI
--------------------

If a package is published to PyPI, you can install directly:

.. code-block:: bash

   pip install rfcnt

Installing from GitHub
----------------------

Install directly from a GitHub release:

.. code-block:: bash

   pip install https://github.com/a-ma72/rainflow/releases/download/rfcnt-0.5.2/rfcnt-0.5.2.tar.gz

Or from the repository:

.. code-block:: bash

   pip install git+https://github.com/a-ma72/rainflow.git

Google Colaboratory
-------------------

To use in Google Colab notebooks:

.. code-block:: python

   !pip install --no-build-isolation --no-deps https://github.com/a-ma72/rainflow/releases/download/rfcnt-0.5.2/rfcnt-0.5.2.tar.gz

   import rfcnt
   rfcnt.tests.examples.example_1()

Running Examples
================

After installing the Python package, run the included examples:

.. code-block:: bash

   python -m rfcnt.run_examples

This will execute all example scripts demonstrating various features of
the package.

Unit Tests
==========

C/C++ Unit Tests
----------------

Build and run the C unit test suite:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_EXPORT_MEX=0
   cmake --build build --target rfc_test --config Release

Then invoke the test executable:

.. code-block:: bash

   # Linux/macOS
   build/test/Release/rfc_test

   # Windows
   build\test\Release\rfc_test.exe

Using CTest
-----------

Alternatively, use CTest to run tests:

.. code-block:: bash

   cd build
   ctest -C Release

Python Unit Tests
-----------------

Run Python tests after installing the package:

.. code-block:: bash

   python -m rfcnt.run_tests

Minimal Build
=============

For embedded systems or microcontrollers, you can build a minimal version
with only core counting functionality:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_MINIMAL=1 -DRFC_EXPORT_PY=0 -DRFC_EXPORT_MEX=0 -DRFC_UNIT_TEST=0
   cmake --build build --config Release

Using COAN for Code Cleanup
----------------------------

Use COAN_ to remove unwanted preprocessor directives and create a clean
minimal version:

.. code-block:: bash

   coan source -DRFC_MINIMAL src/lib/rainflow.c > rainflow_minimal.c

.. _COAN: http://coan2.sourceforge.net/

Custom Feature Selection
=========================

Enable specific features by defining preprocessor macros during the CMake
configuration:

.. code-block:: bash

   cmake -S. -Bbuild \
       -DRFC_TP_SUPPORT=1 \
       -DRFC_HCM_SUPPORT=1 \
       -DRFC_ASTM_SUPPORT=1 \
       -DRFC_DH_SUPPORT=1 \
       -DRFC_DAMAGE_FAST=1

See `features.rst <features.rst>`_ for a complete list of available feature flags.

Platform-Specific Notes
=======================

Windows
-------

- Use Visual Studio 2015 or later
- CMake automatically detects the installed Visual Studio version
- For MinGW, specify: ``-G "MinGW Makefiles"``

Linux
-----

- Install development tools: ``sudo apt-get install build-essential cmake``
- For Python development: ``sudo apt-get install python3-dev``

macOS
-----

- Install Xcode Command Line Tools: ``xcode-select --install``
- Or use Homebrew to install CMake: ``brew install cmake``

Troubleshooting
===============

NumPy Compatibility Issues
--------------------------

If you encounter NumPy API version issues, ensure you're using
``oldest-supported-numpy`` during build:

.. code-block:: bash

   pip install oldest-supported-numpy

CMake Cannot Find MATLAB
-------------------------

If CMake cannot locate your MATLAB installation:

1. Verify MATLAB is installed and accessible
2. Set ``Matlab_ROOT_DIR`` explicitly (see above)
3. Check that your MATLAB version is R2017b or newer

Python Extension Import Errors
-------------------------------

If the Python extension fails to import:

1. Verify NumPy is installed: ``pip install numpy``
2. Check Python architecture matches (32-bit vs 64-bit)
3. Rebuild with ``--no-build-isolation`` flag

See Also
========

- `examples.rst <examples.rst>`_ - Usage after installation
- `features.rst <features.rst>`_ - Available compile-time features
