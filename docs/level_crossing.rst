=======================
Level Crossing Counting
=======================

The level-crossing (LC) histogram counts how often the load series crosses
each class **upper bound** ``u[i] = class_offset + (i+1) * class_width``.
It is independent of rainflow matrix / range-pair residue methods, except
where noted below.

This is **not** ``RFC_RES_RP_DIN45667``, which is a range-pair treatment of
the rainflow residue. See `residue_methods.rst <residue_methods.rst>`_.

Two standards
=============

This library implements two different KGÜZ (Klassengrenzenüberschreitungszählung)
conventions. They disagree as soon as the class grid spans both sides of the
zero-load line.

+------------------+----------------------------------+----------------------------------+
| Aspect           | DIN 45667                        | FVA Merkblatt                    |
+==================+==================================+==================================+
| Counting         | One static slope direction for   | Direction taken from the sign of |
| direction        | **every** class bound            | the class threshold ``u``        |
+------------------+----------------------------------+----------------------------------+
| Sign of signal   | Ignored                          | Ignored (direction comes from    |
|                  |                                  | ``u``, not from the sample)      |
+------------------+----------------------------------+----------------------------------+
| ``u > 0``        | Same global direction            | Positive-going only (away from   |
|                  |                                  | zero into the positive range)    |
+------------------+----------------------------------+----------------------------------+
| ``u < 0``        | Same global direction            | Negative-going only (away from   |
|                  |                                  | zero into the negative range)    |
+------------------+----------------------------------+----------------------------------+
| ``u = 0``        | Same global direction            | Positive branch (rising /        |
|                  |                                  | leaving the zero-load baseline)  |
+------------------+----------------------------------+----------------------------------+
| Context          | Universal statistics of random   | Drivetrain fatigue; load spectra |
|                  | vibrations                       | originating from zero load       |
+------------------+----------------------------------+----------------------------------+
| API              | ``SLOPES_UP`` / ``SLOPES_DOWN``  | ``FVA``. ``DIN45667`` is a       |
|                  | / ``SLOPES_ALL``                 | compatibility alias of ``FVA``   |
+------------------+----------------------------------+----------------------------------+

DIN 45667 — static global direction
-----------------------------------

DIN 45667 is a **universal** statistical classification of random vibrations.
The counting direction is chosen once and applied to **every** class bound,
independent of the algebraic sign of the signal or of ``u``.

Use ``SLOPES_UP`` (positive-going), ``SLOPES_DOWN`` (negative-going), or
``SLOPES_ALL`` (both, stored in the same bin). The library default is
``SLOPES_ALL`` (C ``RFC_FLAGS_COUNT_LC``). Typical DIN practice is one
direction, usually rising (``SLOPES_UP``).

FVA Merkblatt — sign-dependent, from zero load
----------------------------------------------

The FVA drivetrain convention measures exceedances **away from the unloaded
state**. The crossing direction is taken from the sign of the class bound:

* ``u > 0``: count **positive-going** crossings only (leaving zero into the
  positive range).
* ``u < 0``: count **negative-going** crossings only (leaving zero into the
  negative range).
* ``u = 0``: treated as the **positive** branch (rising / leaving the
  zero-load baseline). Counts stay non-negative; the negative branch is not
  sign-flipped.

Use ``LCMethod.FVA`` / ``RFC_LC_COUNT_METHOD_FVA``.
``LCMethod.DIN45667`` / ``RFC_LC_COUNT_METHOD_DIN45667`` is a compatibility
alias of ``FVA`` (historical misnomer: the conversion was previously labelled
DIN 45667).

Slope selection (``lc_method``)
===============================

Python ``LCMethod`` / C ``rfc_lc_count_method`` is an enumeration, not a
flag mask (``0 | 1`` is falling slopes, not both):

+-------+------------------+--------------------------------------------------+
| Value | Name             | Histogram                                        |
+=======+==================+==================================================+
| 0     | ``SLOPES_UP``    | DIN 45667: rising slopes only                    |
+-------+------------------+--------------------------------------------------+
| 1     | ``SLOPES_DOWN``  | DIN 45667: falling slopes only                   |
+-------+------------------+--------------------------------------------------+
| 2     | ``SLOPES_ALL``   | DIN 45667: rising and falling (library default)  |
+-------+------------------+--------------------------------------------------+
| 3     | ``FVA``          | FVA: both internally; convert on ``RFC_lc_get``  |
+-------+------------------+--------------------------------------------------+

Default ``2`` matches C ``RFC_FLAGS_COUNT_LC``. Combined counting (``2``)
does **not** match FVA, which assigns the direction from the **sign of the
class boundary**.

FVA conversion
==============

The library still increments ``ctx.lc`` with **both** directions
(``n_ges``). ``RFC_lc_get`` / Python ``lc`` then convert.

For each level ``u``::

   delta(u) = 1[x_end > u] - 1[x_start > u]   in {-1, 0, +1}

``delta`` is nonzero only for ``min(x_start, x_end) < u < max(...)``
(strict inequalities: levels exactly equal to an endpoint are unaffected)::

   n_plus(u)  = (n_ges(u) + delta(u)) / 2
   n_minus(u) = (n_ges(u) - delta(u)) / 2
   n_FVA(u)   = n_plus(u)   if u >= 0
   n_FVA(u)   = n_minus(u)  if u <  0

This is equivalent to counting ``SLOPES_UP`` on non-negative bounds and
``SLOPES_DOWN`` on negative bounds. Division by two must be exact. A
remainder is ``RFC_ERROR_DATA_INCONSISTENT`` (Python: runtime error);
values are **not** rounded.

``x_start`` / ``x_end`` are the first and last **fed samples**, not turning
points. Residue finalizers do not add further LC increments (already
stripped for all LC methods). With ``FVA``, ``RFC_lc_get`` includes the
interim turning-point slope (if the context is still ``BUSY_INTERIM``) so
live ``RFC.lc`` matches ``RFC.lc_as(...)`` and the histogram after
``finalize()``.

Prerequisite
------------

Use ``enforce_margin=True`` (Python default) so the first and last samples
are turning points. Otherwise hysteresis can make ``n_ges`` disagree with
the raw endpoints and check P2 fails.

C API: set ``ctx.lc_count_method = RFC_LC_COUNT_METHOD_FVA`` after
``RFC_init`` (keep ``RFC_FLAGS_COUNT_LC``), or call
``RFC_lc_convert_fva`` on an existing combined histogram.
``RFC_lc_convert_din45667`` is a compatibility wrapper for
``RFC_lc_convert_fva``.

``RFC_lc_from_rfm`` and ``RFC_lc_from_residue`` reconstruct a histogram from
closed cycles or a residue using ``RFC_FLAGS_COUNT_LC_UP`` / ``_DN`` only
(DIN static direction). They do **not** apply FVA conversion, which needs
the series endpoints ``x_start`` / ``x_end``. Use ``RFC_lc_get`` after
feeding the time series for FVA results.

MATLAB MEX ``rfc()`` returns the combined DIN histogram (both slopes), the
same as C default ``RFC_FLAGS_COUNT_LC``. It has no ``lc_method`` argument.

Checks after conversion
-----------------------

P1
   ``n_FVA[i]`` is an integer and ``>= 0``.

P2
   ``n_ges[i]`` is odd if and only if ``lo < u[i] < hi``. An odd count
   outside that interval usually means the counter treated a sample exactly
   on a class boundary inconsistently, or missed a first/last turning point.

P3
   From ``u = 0``, ``n_FVA`` is monotone toward each end (non-decreasing
   toward zero on the negative side, non-increasing away from zero on the
   positive side). An increase going outward is always an error.

P4
   The maximum is at the bin(s) adjacent to ``u = 0`` (no class needs
   ``u == 0`` exactly; the test vectors peak at ±0.5).

Test vectors
------------

Class upper bounds on half-integers.

**T1** ``x = [-3, +2, -1]`` (``x_end > x_start``; affected levels ``-2.5``,
``-1.5``, both negative => ``delta = +1``):

::

   u      = [-2.5, -1.5, -0.5, +0.5, +1.5]
   n_ges  = [   1,    1,    2,    2,    2]
   n_FVA  = [   0,    0,    1,    1,    1]

DIN ``SLOPES_UP`` on the same series is ``[1, 1, 1, 1, 1]``.

**T2** ``x = [+1, -2, +3]`` (affected ``+1.5``, ``+2.5``, both positive =>
``delta = +1``):

::

   u      = [-1.5, -0.5, +0.5, +1.5, +2.5]
   n_ges  = [   2,    2,    2,    1,    1]
   n_FVA  = [   1,    1,    1,    1,    1]

Conservative alternative (not implemented)
------------------------------------------

Appending ``x_start`` at the end of the series before counting makes
``delta`` identically zero. Then ``n_ges`` is even on every level and
``n_ges / 2`` is already the FVA result. The extra limb raises the count
by 1 on the affected levels (conservative). Do not use rounding as a
substitute for the ``delta`` case distinction.

See also
========

- `features.rst <features.rst>`_ - Histogram overview
- `residue_methods.rst <residue_methods.rst>`_ - ``RFC_RES_RP_DIN45667``
- `algorithm.rst <algorithm.rst>`_ - 4-point counting
- `references.rst <references.rst>`_ - DIN 45667 and FVA-Richtlinie citations
