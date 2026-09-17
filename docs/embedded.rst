=======================
Embedded Rainflow Engine
=======================

Overview
========

``src/embedded/`` is a **separate** fixed-point rainflow engine for small
devices (MCUs and boards without an FPU). Its public API uses the ``RF_*``
prefix (``RF_Init``, ``RF_ProcessSample``, ``RF_FlushResiduumRepeated``, …).

It is **not** ``RFC_MINIMAL``.

+---------------------------+--------------------------------------------+--------------------------------------------+
|                           | ``src/lib`` + ``RFC_MINIMAL``              | ``src/embedded`` (``RF_*``)                |
+===========================+============================================+============================================+
| What it is                | Compile-time slice of the full ``RFC_*``   | Independent port with its own API          |
|                           | library (``unifdef`` / feature flags)      |                                            |
+---------------------------+--------------------------------------------+--------------------------------------------+
| Sample type               | ``double`` (``rfc_value_t``)               | Q28.4 fixed-point (``RF_Value_t``)         |
+---------------------------+--------------------------------------------+--------------------------------------------+
| Heap                      | Allocator / optional buffers               | None — context is fully static             |
+---------------------------+--------------------------------------------+--------------------------------------------+
| Working storage           | Residue, RFM, optional TP / DH             | Residue stack only                         |
+---------------------------+--------------------------------------------+--------------------------------------------+
| Counting path             | Floating-point classes and damage          | Integer class indices + compile-time LUT   |
+---------------------------+--------------------------------------------+--------------------------------------------+
| Residue method            | Several ``RFC_RES_*`` options              | Repeated residue (commit or predict)       |
+---------------------------+--------------------------------------------+--------------------------------------------+

The 4-point closure rule is the same as in ``src/lib`` (DIN 45667 / FVA).
The implementations are independent: this engine is conceptually aligned
with that method, not a compile-flag extract of ``rainflow.c``.

Include the headers from ``src/embedded/`` (or as ``embedded/rainflow.h``
via ``src/``) so they never collide with ``src/lib/rainflow.h``.

When to use this engine
=======================

Use ``RF_*`` when:

- The target has **no FPU** (or software float is too expensive)
- RAM must stay bounded and **static** (no ``malloc``)
- The input is a **never-ending stream** (no full time-series buffer)
- You only need **accumulated damage** and an optional range-pair histogram

Use ``RFC_MINIMAL`` (or the full ``RFC_*`` library) when you need a rainflow
matrix, turning points, HCM/ASTM, damage history, or floating-point classes.
See `minimal_build.rst <minimal_build.rst>`_.

Two-stage processing
====================

Each incoming sample is processed in two stages:

1. **Q28.4 hysteresis / turning-point filter** on the raw sample.
   A turning point is confirmed only after it exceeds the configured
   hysteresis. Noise is suppressed **before** quantization.

2. **Integer classification.** A confirmed turning point is mapped to a
   class index (``0 .. RF_NUM_CLASSES-1``) and pushed onto the residue
   stack. The 4-point algorithm then uses **only integers** — no
   fixed-point arithmetic on the hot counting path.

Closed cycles do not emit from/to pairs. Amplitude → damage is a table
lookup (``RF_DAMAGE_LUT``) at cycle closure.

Hysteresis precondition
-----------------------

``RF_Init()`` requires ``hysteresis >= RF_CLASS_WIDTH`` (same Q28.4 units).
Only then are two consecutive confirmed turning points guaranteed to land
in different classes. A smaller hysteresis can break residue alternation
and the correctness argument of ``RF_FlushResiduumRepeated()``. Init
returns ``RF_ERR_INVALID_CONFIG`` if the precondition fails.

Residue and memory
==================

Only the classified residue stack is kept. Maximum depth is
``RF_MAX_RESIDUUM = 2 * RF_NUM_CLASSES`` (each class can appear at most
twice as an open flank). That bound holds only when the hysteresis
precondition is met.

There is no time-series buffer, no turning-point list, and no rainflow
matrix. The context (``RF_Ctx_t``) is fully static and **not thread-safe**.

Public API
==========

.. code-block:: c

   RF_Status_t   RF_Init(RF_Ctx_t *ctx, RF_Value_t hysteresis);
   RF_Status_t   RF_ProcessSample(RF_Ctx_t *ctx, RF_Value_t sample,
                                  uint32_t *out_damage_increment);
   RF_Status_t   RF_FlushResiduumRepeated(RF_Ctx_t *ctx, bool commit,
                                          uint32_t *out_damage_increment);
   uint16_t      RF_GetResiduumCount(const RF_Ctx_t *ctx);
   RF_Damage96_t RF_GetAccumulatedDamage(const RF_Ctx_t *ctx);
   double        RF_GetAccumulatedDamageDouble(const RF_Ctx_t *ctx);
   const uint32_t *RF_GetRangePairCounts(const RF_Ctx_t *ctx);

``RF_ProcessSample``
   Feed one Q28.4 sample. If a turning point is confirmed and cycles
   close, their LUT damage is added to the 96-bit accumulator and to
   ``*out_damage_increment`` (may be ``NULL``). Range-pair counts
   (``rp_counts[range]``) increment for each closed cycle.

``RF_FlushResiduumRepeated``
   Apply the same *repeated residue* idea as ``RFC_RES_REPEATED`` in
   ``src/lib``: notionally append the residue to itself and extract
   extra full cycles.

   - ``commit == false`` (**Predict**): compute the increment, leave the
     context unchanged (including residue and accumulated damage).
   - ``commit == true``: apply the increment, reduce the residue to the
     last point, and continue the stream from there.

   Because the stream never really ends, call this periodically with
   ``commit == true`` (for example every N seconds) so open remainder
   cycles are not left uncounted forever.

``RF_GetAccumulatedDamage``
   Exact 96-bit total of **applied** increments (sample closures plus
   committed flushes). Predict calls do not contribute.

``RF_GetAccumulatedDamageDouble``
   Convenience for host tools. **Do not** call it on targets without an
   FPU (it pulls in software float). Above ``2^53`` a ``double`` cannot
   represent every integer; use the 96-bit getter for the exact total.

``RF_GetRangePairCounts``
   Histogram over class difference (same index as ``RF_DAMAGE_LUT``).
   Optional in production: the ``rp_counts`` field may be omitted to
   save RAM if the histogram is not needed.

Example
=======

.. code-block:: c

   #include "rainflow.h"   /* from src/embedded/, not src/lib/ */

   RF_Ctx_t ctx;
   uint32_t increment;

   if (RF_Init(&ctx, RF_CLASS_WIDTH_FIXED) != RF_OK)
   {
       /* hysteresis must be >= RF_CLASS_WIDTH */
       return;
   }

   /* samples already in Q28.4, or convert with RF_DOUBLE_TO_FIXED */
   RF_ProcessSample(&ctx, sample, &increment);

   /* Preview residue damage without changing state */
   RF_FlushResiduumRepeated(&ctx, false, &increment);

   /* Periodically commit remainder cycles */
   RF_FlushResiduumRepeated(&ctx, true, NULL);

   RF_Damage96_t D = RF_GetAccumulatedDamage(&ctx);
   const uint32_t *rp = RF_GetRangePairCounts(&ctx);

Configuration and damage LUT
============================

Application-specific class grid and the damage table live in
``rainflow_config.h`` / ``rainflow_config.c``:

- ``RF_NUM_CLASSES``, ``RF_CLASS_WIDTH``, ``RF_CLASS_MIN``
- ``RF_FIXED_SHIFT`` (Q28.4: 4 fractional bits)
- ``RF_DAMAGE_LUT[i]`` — damage of one closed cycle with class
  difference ``i``

Regenerate the example LUT from a Wöhler exponent ``k``::

   python src/embedded/generate_rainflow_config.py --num-classes 128 --k 5 \
       --output src/embedded/rainflow_config.c

``RF_DAMAGE_LUT[i]`` is ``round(i ** k)`` (so index 0 is always 0).
Replace this with your material S-N curve for production firmware.

Host tests
==========

With the main CMake tree (``RFC_UNIT_TEST=ON``, the default for a standalone
build)::

   cmake --build build --target rf_embedded_test rf_embedded_flush_test
   ctest -C Release -R rf_embedded -V

Or from ``src/embedded/test`` (assert-based, no Greatest)::

   gcc -std=c11 -Wall -Wextra -Wpedantic -I.. ../rainflow.c test_rainflow.c -o test_rainflow
   ./test_rainflow

   gcc -O2 -Wall -Wextra -o test_RF_FlushResiduumRepeated test_RF_FlushResiduumRepeated.c
   ./test_RF_FlushResiduumRepeated

``test_rainflow.c`` supplies its own identity ``RF_DAMAGE_LUT`` and must
**not** be linked with ``rainflow_config.c``.

See Also
========

- `minimal_build.rst <minimal_build.rst>`_ - ``RFC_MINIMAL`` slice of the full library
- `algorithm.rst <algorithm.rst>`_ - 4-point method (same closure rule)
- `residue_methods.rst <residue_methods.rst>`_ - ``RFC_RES_REPEATED`` in ``src/lib``
- `features.rst <features.rst>`_ - Compile-time flags for ``RFC_*``

Source files: ``src/embedded/rainflow.h``, ``src/embedded/rainflow.c``,
``src/embedded/rainflow_config.h``.
