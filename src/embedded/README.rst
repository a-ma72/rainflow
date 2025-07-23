Embedded rainflow engine (``RF_*``)
===================================

This directory is a **separate** fixed-point rainflow engine for small
devices (MCUs, boards without an FPU). It is **not** ``RFC_MINIMAL``.

Full documentation lives in the main docs:

`docs/embedded.rst <../../docs/embedded.rst>`_

``src/lib/minimal/`` is an ``unifdef`` extract of the full ``RFC_*``
library. This port has its own API (``RF_Init``, ``RF_ProcessSample``,
``RF_FlushResiduumRepeated``, …). Include the headers from this
directory so they never collide with ``src/lib/rainflow.h``.

Layout
------

* ``rainflow.c`` / ``rainflow.h`` — engine
* ``rainflow_config.h`` / ``rainflow_config.c`` — example class grid and LUT
* ``rainflow_damage_lut_meta.h`` — Wöhler exponent used for the LUT
* ``generate_rainflow_config.py`` — regenerates ``rainflow_config.c``
* ``test/`` — host tests (assert-based, no Greatest; CMake targets
  ``rf_embedded_test`` and ``rf_embedded_flush_test``)
