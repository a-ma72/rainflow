=============================================
Rainflow Counting Algorithm (C99 compliant)
=============================================

.. image:: https://github.com/a-ma72/rainflow/actions/workflows/run_test.yml/badge.svg
   :target: https://github.com/a-ma72/rainflow/actions/workflows/run_test.yml
   :alt: tests

A robust, modular implementation of the Rainflow Counting Algorithm using
the 4-point method for fatigue analysis. This library is C99 compliant
and provides bindings for Python and MATLAB.

Quick Start
===========

Python
------

Install from PyPI:

.. code-block:: bash

   pip install rfcnt

Basic usage:

.. code-block:: python

   import numpy as np
   import rfcnt

   # Your measurement data
   data = np.array([0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0])

   # Perform rainflow counting
   result = rfcnt.rfc(data, class_width=0.5)

   print(f"Damage: {result['damage']}")
   print(f"Cycles: {result['rp']}")

C/C++
-----

.. code-block:: c

   #include "rainflow.h"

   rfc_ctx_s ctx;
   double data[] = {0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0};

   RFC_init(&ctx, 100, 0.5, 0.0, 0.5, RFC_FLAGS_DEFAULT);
   RFC_feed(&ctx, data, 7);
   RFC_finalize(&ctx, RFC_RES_REPEATED);

   double damage;
   RFC_damage(&ctx, &damage);
   RFC_deinit(&ctx);

Key Features
============

- **Modular architecture** - Select features at compile time
- **Multiple counting methods** - 4-point, HCM (Clormann/Seeger), ASTM
- **Streaming capability** - Process data sample-by-sample or in chunks
- **Fatigue analysis** - Wöhler curves, Miner's rule (4 variants), damage history
- **Flexible histograms** - Rainflow matrix, level crossing, range pairs
- **Language bindings** - Python (NumPy), MATLAB, C/C++
- **Standards compliant** - ASTM E 1049, DIN 45667, FKM guidelines

Documentation
=============

.. contents:: Table of Contents
   :depth: 2
   :local:

Detailed Documentation
----------------------

Explore the full documentation in the ``docs/`` folder:

`Documentation Index <docs/index.rst>`_
   Complete documentation index with all topics organized by category

Core Documentation:

`Algorithm <docs/algorithm.rst>`_
   Detailed explanation of the rainflow counting algorithm, including the
   4-point method, hysteresis filtering, and residue handling.

`ASTM Method <docs/astm_method.rst>`_
   ASTM E 1049 counting method vs. 4-point algorithm - why it's more efficient

`Features <docs/features.rst>`_
   Complete feature list with compile-time options, counting methods,
   fatigue analysis capabilities, and language bindings.

`Installation <docs/installation.rst>`_
   Build instructions for all platforms (Windows, Linux, macOS) with
   Python, MATLAB, and C/C++ configurations.

`Examples <docs/examples.rst>`_
   Practical code examples in Python, C/C++, and MATLAB demonstrating
   various use cases and workflows.

Advanced Topics:

`C++ Wrapper <docs/cpp_wrapper.rst>`_
   Modern C++ interface with RAII, templates, and STL integration

`Damage History <docs/damage_history.rst>`_
   Track cumulative damage over time for predictive maintenance

`Delegates <docs/delegates.rst>`_
   Custom behavior via function pointers and callbacks

`Turning Points <docs/turning_points.rst>`_
   External storage and management of turning points

`Residue Methods <docs/residue_methods.rst>`_
   Different approaches for handling unclosed cycles

`Spread Damage <docs/spread_damage.rst>`_
   Transient damage distribution methods explained

`Minimal Build <docs/minimal_build.rst>`_
   RFC_MINIMAL for embedded systems and microcontrollers

`TP Prune <docs/tp_prune.rst>`_
   RFC_tp_prune() logic for memory management

`References <docs/references.rst>`_
   Bibliography and citations for standards, publications, and resources.

Quick Links
-----------

- `GitHub Repository <https://github.com/a-ma72/rainflow>`_
- `Issue Tracker <https://github.com/a-ma72/rainflow/issues>`_
- `Releases <https://github.com/a-ma72/rainflow/releases>`_

What is Rainflow Counting?
===========================

Rainflow Counting is a standardized method for analyzing fatigue in materials
subjected to variable amplitude loading. The algorithm extracts discrete
load cycles from complex time-varying stress-strain histories.

The Four-Point Method
----------------------

The core algorithm examines sequences of four consecutive turning points
(A, B, C, D) to identify closed cycles::

                     * D
                    / \
             B *<--/          A cycle B-C is closed if:
              / \ /           min(B,C) >= min(A,D) &&
             /   * C          max(B,C) <= max(A,D)
          \ /
           * A

When closed, the cycle B-C is:

1. Counted in the histogram
2. Removed from the residue
3. Assigned damage based on the configured Wöhler curve

See `docs/algorithm.rst <docs/algorithm.rst>`_ for detailed explanation.

Installation
============

Build from Sources
------------------

Using CMake (all platforms):

.. code-block:: bash

   cmake -S. -Bbuild -G "Visual Studio 16 2019"
   cmake --build build --config Release

Or on Linux/macOS:

.. code-block:: bash

   cmake -S. -Bbuild
   make -C build

Python Wheel
------------

Build a Python wheel package:

.. code-block:: bash

   cd src/python
   pip install setuptools build wheel oldest-supported-numpy
   python -m build -nw

For detailed build instructions, including MATLAB integration and custom
feature selection, see `docs/installation.rst <docs/installation.rst>`_.

Testing
=======

Run the unit tests to verify your installation:

C/C++ Tests
-----------

.. code-block:: bash

   cmake --build build --target rfc_test --config Release
   build/test/Release/rfc_test      # Linux/macOS
   build\test\Release\rfc_test.exe  # Windows

Or use CTest:

.. code-block:: bash

   cd build
   ctest -C Release

Python Tests
------------

.. code-block:: bash

   python -m rfcnt.run_tests

Run Examples
------------

.. code-block:: bash

   python -m rfcnt.run_examples

Current Status
==============

- **Version**: 0.5.2
- **C Standard**: C99 compliant
- **Test Status**: |tests|
- **Languages**: C, C++, Python, MATLAB
- **Platforms**: Windows, Linux, macOS

.. |tests| image:: https://github.com/a-ma72/rainflow/actions/workflows/run_test.yml/badge.svg
   :target: https://github.com/a-ma72/rainflow/actions/workflows/run_test.yml

Applications
============

This library is suitable for:

- **Fatigue life prediction** in mechanical engineering
- **Structural health monitoring** of bridges, vehicles, aircraft
- **Materials testing** and characterization
- **Load spectrum analysis** for component design
- **Real-time monitoring** in embedded systems

Standards Compliance
====================

The implementation follows these standards:

- **ASTM E 1049** (2011) - Standard Practices for Cycle Counting
- **DIN 45667-1:2025-12** - Fatigue analysis of structures - Part 1: Rainflow counting method
- **DIN 45667-2:2025-12** - Fatigue analysis of structures - Part 2: Damage accumulation
- **FKM Guidelines** - Analytical strength assessment
- **FVA-Richtlinie** - Counting methods for drive technology

See `docs/references.rst <docs/references.rst>`_ for full bibliography.

Architecture
============

Modular Design
--------------

The package uses a two-layer architecture:

**Core C Library** (``src/lib/rainflow.c``)
   Pure C99 implementation with optional feature flags for customization.
   Suitable for embedded systems and microcontrollers when built with
   ``RFC_MINIMAL``.

**C++ Wrapper** (``src/lib/rainflow.hpp``)
   Modern C++ interface with templates, containers, and RAII patterns.
   Provides object-oriented access and easier integration with C++ projects.

**Language Bindings**
   - Python extension (``src/python/src/rfcnt.cpp``) with NumPy support
   - MATLAB MEX interface for use in MATLAB/Simulink

Feature Flags
-------------

Customize the build with compile-time flags:

- ``RFC_MINIMAL`` - Core counting only
- ``RFC_TP_SUPPORT`` - Turning point storage
- ``RFC_HCM_SUPPORT`` - HCM algorithm
- ``RFC_ASTM_SUPPORT`` - ASTM variant
- ``RFC_DH_SUPPORT`` - Damage history
- ``RFC_DAMAGE_FAST`` - Lookup tables for speed

See `docs/features.rst <docs/features.rst>`_ for complete list.

Contributing
============

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

Report issues at: https://github.com/a-ma72/rainflow/issues

License
=======

See the LICENSE file in the repository for licensing information.

Citation
========

If you use this software in your research, please cite:

.. code-block:: bibtex

   @software{rainflow,
     author = {Andreas Martin},
     title = {Rainflow Counting Algorithm},
     year = {2026},
     url = {https://github.com/a-ma72/rainflow},
     version = {0.5.2}
   }

Acknowledgments
===============

This implementation is based on the work described in the references
section. Special thanks to the authors of ASTM E 1049 and the HCM method
(Clormann/Seeger).

Contact
=======

For questions, suggestions, or collaboration:

- GitHub Issues: https://github.com/a-ma72/rainflow/issues
- Repository: https://github.com/a-ma72/rainflow

---

**Documentation Structure**:

.. code-block:: text

   rainflow/
   ├── README.rst                   ← Readme
   └── docs/
       ├── index.rst                ← Documentation index (start here)
       ├── algorithm.rst            ← Algorithm explanation
       ├── features.rst             ← Feature descriptions
       ├── installation.rst         ← Build instructions
       ├── examples.rst             ← Code examples
       ├── references.rst           ← Bibliography
       ├── astm_method.rst          ← ASTM E1049 vs 4-point method
       ├── cpp_wrapper.rst          ← C++ wrapper guide
       ├── damage_history.rst       ← Damage tracking over time
       ├── delegates.rst            ← Custom callbacks & extensibility
       ├── turning_points.rst       ← Turning point storage & management
       ├── residue_methods.rst      ← Residue handling methods
       ├── spread_damage.rst        ← Damage spreading methods
       ├── minimal_build.rst        ← Embedded systems build
       └── tp_prune.rst             ← Memory management with pruning

Start exploring with `docs/index.rst <docs/index.rst>`_ for organized navigation,
`docs/algorithm.rst <docs/algorithm.rst>`_ for algorithm details, or jump
directly to `docs/examples.rst <docs/examples.rst>`_ for practical usage.
