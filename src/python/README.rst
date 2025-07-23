Rainflow Python Package (rfcnt)
===============================

Summary
-------

The **rfcnt** package provides a high-performance, standards-compliant implementation of the rainflow counting algorithm for fatigue analysis, with bindings for Python, C/C++, and MATLAB. It supports the 4-point method (DIN 45667), the ASTM E 1049 3-point method, HCM (Clormann/Seeger), and advanced residue handling. This README summarizes the package for users and developers. For full documentation, see the `/docs` folder.

Features
--------

- Fast C core with Python and MATLAB wrappers
- 4-point, HCM, and ASTM counting methods
- Flexible residue processing (DIN, ASTM, repeated, ignore, HCM)
- Streaming and batch processing
- Dynamic class management (auto-resize, offset)
- Compile-time feature selection (minimal, delegates, damage history, etc.)
- Amplitude transformation (FKM Haigh / ``rfcnt.at_transform``)
- Unit tests and real-world examples

Algorithm Overview
------------------

Rainflow counting extracts closed cycles from load histories in four main steps:

1. **Hysteresis Filtering**: Removes small oscillations below a threshold
2. **Peak-Valley Filtering**: Identifies turning points
3. **Discretization**: Maps values to classes (bins)
4. **4-Point Counting**: Detects closed cycles using the 4-point pattern

A cycle B-C is closed if:

    min(B, C) >= min(A, D) and max(B, C) <= max(A, D)

Closed cycles are counted and removed; residue is handled per user choice.

Installation
------------

**Prerequisites:**
- Python 3.6+
- NumPy (>=2.0)
- CMake 3.16+
- C/C++ compiler (GCC, Clang, MSVC)

**Build and Install:**

    cmake -S. -Bbuild -DRFC_EXPORT_PY=1
    cmake --build build --target build_wheel_isolated --config Release
    pip install dist/rfcnt-*.whl

Or install from source:

    pip install . --no-build-isolation --no-deps

**MATLAB:**
See docs/installation.rst for MATLAB integration.

Quick Start
-----------

    import numpy as np
    import rfcnt

    data = np.array([0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0])
    result = rfcnt.rfc(data, class_width=0.5, class_count=10)
    print("Total damage:", result['damage'])
    print("Range pairs:\n", result['rp'])

Stateful counting (chunked feed)
--------------------------------

``rfc()`` is one-shot. For a persistent counter, construct ``RFC``, call
``feed()`` one or more times, then ``finalize()``:

    from rfcnt import RFC, ResidualMethod

    rf = RFC(class_width=0.5, class_count=100, wl={"sx": 1e3, "nx": 1e7, "k": 5})
    rf.feed(chunk1)
    rf.feed(chunk2)
    print(rf.damage_as(ResidualMethod.REPEATED))  # preview; object stays open
    rf.finalize(ResidualMethod.REPEATED)
    print(rf.damage, rf.residue, rf.res_raw)
    print(rf.tp)  # (n, 4): pos, value, damage, adj_pos

Damage history is not supported on ``RFC`` (``spread_damage`` must stay
``SDMethod.NONE``); use one-shot ``rfc()`` when you need a damage history. After
``finalize()``, further ``feed()`` calls raise. Start a new ``RFC`` for a new series.
``damage`` is live closed-cycle damage; ``damage_as(method)`` is that value plus
residue processed as if ``finalize(method)`` had been called. The same split
applies to ``rp`` / ``lc`` / ``rfm`` and ``rp_as`` / ``lc_as`` / ``rfm_as``.
``tp`` is live and may grow when ``finalize()`` promotes the interim point.
``res_raw`` is the open residue (4-point closed-cycle strip, same as
``rfc()["res_raw"]``): live until ``finalize()``, then an isolated read-only
snapshot. ``wl_miner_consistent`` is live (impaired Wöhler parameters from
closed cycles so far), like ``damage``, not a ``damage_as`` preview. After
``finalize()`` it matches ``rfc()["wl_miner_consistent"]``.

Advanced Usage
--------------

- Use Wöhler curve parameters for fatigue life prediction
- Select counting method: `use_HCM`, `use_ASTM`
- Select level-crossing slopes: `lc_method` (default both, like C ``RFC_FLAGS_COUNT_LC``).
  ``LCMethod.SLOPES_UP`` / ``SLOPES_DOWN`` / ``SLOPES_ALL`` are DIN 45667
  (static global direction). ``LCMethod.FVA`` converts the combined histogram
  to the FVA Merkblatt convention (sign-dependent crossings from a zero-load
  baseline); residue methods do not change that ``lc``.
  ``LCMethod.DIN45667`` is a compatibility alias of ``FVA``.
- Control residue handling: `residual_method` (``rfc()``) or ``finalize()`` / ``damage_as()`` (``RFC``)
- Enable damage history on one-shot ``rfc()``: `spread_damage` (not supported on ``RFC``)
- Integrate with real-time or embedded systems (RFC_MINIMAL)

Example:

    result = rfcnt.rfc(
        data,
        class_width=5.0,
        use_HCM=True,
        spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
        hysteresis=2.5
    )

Documentation
-------------

- Algorithm details: docs/algorithm.rst
- Features: docs/features.rst
- Examples: docs/examples.rst
- Installation: docs/installation.rst
- References: docs/references.rst

References
----------

[1] ASTM E 1049, "Standard Practices for Cycle Counting in Fatigue Analysis", ASTM International, 2011.

[2] U.H. Clormann, T. Seeger, "Rainflow - HCM / Ein Hysteresisschleifen-Zaehlalgorithmus...", TU Darmstadt, 1985.

[3] FVA-Richtlinie, 2010. https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf

[4] Siemens PLM, "Rainflow Counting", 2018. https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093

For a full bibliography, see docs/references.rst.

License
-------

See LICENSE file.

Contact
-------

For questions, bug reports, or contributions, see the project repository or contact the maintainer.

