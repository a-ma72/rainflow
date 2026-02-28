====================================
Rainflow Python Package (rfcnt)
====================================

Summary
=======

The **rfcnt** package provides a high-performance, standards-compliant implementation of the rainflow counting algorithm for fatigue analysis, with bindings for Python, C/C++, and MATLAB. It supports the 4-point method (ASTM E 1049, DIN 45667), HCM (Clormann/Seeger), and advanced residue handling. This README summarizes the package for users and developers. For full documentation, see the `/docs` folder.

Features
========
- Fast C core with Python and MATLAB wrappers
- 4-point, HCM, and ASTM counting methods
- Flexible residue processing (DIN, ASTM, repeated, ignore, HCM)
- Streaming and batch processing
- Dynamic class management (auto-resize, offset)
- Compile-time feature selection (minimal, delegates, damage history, etc.)
- Unit tests and real-world examples

Algorithm Overview
==================
Rainflow counting extracts closed cycles from load histories in four main steps:

1. **Hysteresis Filtering**: Removes small oscillations below a threshold
2. **Peak-Valley Filtering**: Identifies turning points
3. **Discretization**: Maps values to classes (bins)
4. **4-Point Counting**: Detects closed cycles using the 4-point pattern

A cycle B-C is closed if:

    min(B, C) >= min(A, D) and max(B, C) <= max(A, D)

Closed cycles are counted and removed; residue is handled per user choice.

Installation
============

**Prerequisites:**
- Python 3.6+
- NumPy (>=2.0)
- CMake 3.10+
- C/C++ compiler (GCC, Clang, MSVC)

**Build and Install:**

.. code-block:: bash

    cmake -S. -Bbuild -DRFC_EXPORT_PY=1
    cmake --build build --target build_wheel_isolated --config Release
    pip install dist/rfcnt-*.whl

Or install from source:

.. code-block:: bash

    pip install . --no-build-isolation --no-deps

**MATLAB:**
See docs/installation.rst for MATLAB integration.

Quick Start
===========

.. code-block:: python

    import numpy as np
    import rfcnt

    data = np.array([0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0])
    result = rfcnt.rfc(data, class_width=0.5, class_count=10)
    print("Total damage:", result['damage'])
    print("Range pairs:\n", result['rp'])

Advanced Usage
==============

- Use Wöhler curve parameters for fatigue life prediction
- Select counting method: `use_HCM`, `use_ASTM`
- Control residue handling: `residual_method`
- Enable damage history: `spread_damage`
- Integrate with real-time or embedded systems (RFC_MINIMAL)

Example:

.. code-block:: python

    result = rfcnt.rfc(
        data,
        class_width=5.0,
        use_HCM=True,
        spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
        hysteresis=2.5
    )

Documentation
=============

- Algorithm details: docs/algorithm.rst
- Features: docs/features.rst
- Examples: docs/examples.rst
- Installation: docs/installation.rst
- References: docs/references.rst

References
==========

[1] ASTM E 1049, "Standard Practices for Cycle Counting in Fatigue Analysis", ASTM International, 2011.
[2] U.H. Clormann, T. Seeger, "Rainflow - HCM / Ein Hysteresisschleifen-Zaehlalgorithmus...", TU Darmstadt, 1985.
[3] FVA-Richtlinie, 2010. https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf
[4] Siemens PLM, "Rainflow Counting", 2018. https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093

For a full bibliography, see docs/references.rst.

License
=======
See LICENSE file.

Contact
=======
For questions, bug reports, or contributions, see the project repository or contact the maintainer.
