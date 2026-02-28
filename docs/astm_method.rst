=============================================
ASTM E 1049 vs. 4-Point Counting Method
=============================================

Overview
========

This library implements two rainflow counting algorithms:

- **4-Point Method** (DIN 45667 / FVA-Richtlinie) - Examines 4 consecutive turning points
- **ASTM 3-Point Method** (ASTM E 1049) - Examines 3 consecutive turning points

Both methods produce equivalent results for complete load histories, but they differ
in their approach to cycle detection and handling of the first turning point.

ASTM E 1049 Standard
=====================

Background
----------

ASTM E 1049 - "Standard Practices for Cycle Counting in Fatigue Analysis"
(2011 revision) is a comprehensive standard that describes multiple cycle
counting methods:

- Range counting
- Peak counting
- Level crossing counting
- **Rainflow counting** (3-point method)

The rainflow counting method in ASTM E 1049 uses a 3-point algorithm with
special handling for cycles that include the first turning point.

4-Point Method (DIN 45667)
===========================

How It Works
------------

The 4-point algorithm examines four consecutive turning points (A, B, C, D)
to determine if points B and C form a closed cycle::

                     * D
                    / \
             B *<--/          Closed if:
              / \ /           min(B,C) >= min(A,D) &&
             /   * C          max(B,C) <= max(A,D)
          \ /
           * A

**Implementation** (from ``cycle_find_4ptm()``):

.. code-block:: c

   while( rfc_ctx->residue_cnt >= 4 )
   {
       size_t idx = rfc_ctx->residue_cnt - 4;

       unsigned A = rfc_ctx->residue[idx+0].cls;
       unsigned B = rfc_ctx->residue[idx+1].cls;
       unsigned C = rfc_ctx->residue[idx+2].cls;
       unsigned D = rfc_ctx->residue[idx+3].cls;

       // Sort B,C and A,D for min/max comparison
       if( B > C ) { swap(B, C); }
       if( A > D ) { swap(A, D); }

       // Check closure condition
       if( A <= B && C <= D )
       {
           // Cycle B-C is closed
           // Count cycle and remove B,C from residue
           ...
       }
       else break;
   }

**Algorithm Steps:**

1. Read the last four turning points from residue: A, B, C, D
2. Sort to get min/max: min(B,C), max(B,C), min(A,D), max(A,D)
3. Check closure condition: ``min(A,D) <= min(B,C) && max(B,C) <= max(A,D)``
4. If closed:

   - Count cycle B-C (full cycle)
   - Remove B and C from residue
   - Repeat from step 1

5. If not closed, stop checking (new turning point needed)

**Characteristics:**

- Always counts closed cycles as **full cycles**
- Simple closure condition
- No special handling for first turning point

ASTM 3-Point Method
====================

How It Works
------------

The ASTM method examines only three consecutive turning points (A, B, C)
and compares their ranges::


             B *              Range Y = |A - B|
              / \             Range X = |B - C|
             /   \
          \ /     \ /         Closed if: X >= Y
           * A     * C

**Implementation** (from ``cycle_find_astm()``):

.. code-block:: c

   while( rfc_ctx->residue_cnt >= 3 )
   {
       size_t idx = rfc_ctx->residue_cnt - 3;

       unsigned A = rfc_ctx->residue[idx+0].cls;
       unsigned B = rfc_ctx->residue[idx+1].cls;
       unsigned C = rfc_ctx->residue[idx+2].cls;
       unsigned Y = abs( (int)A - (int)B );  // Range Y
       unsigned X = abs( (int)B - (int)C );  // Range X

       // Check closure condition
       if( X >= Y )
       {
           rfc_value_tuple_s *Z = &rfc_ctx->residue[0];  // First TP

           // Does range Y include the first turning point Z?
           if( (Z->cls >= A && Z->cls <= B) ||
               (Z->cls >= B && Z->cls <= A) )
           {
               // Count as HALF cycle
               cycle_process_counts(...);
               // Remove only A (one point)
               residue_cnt--;
           }
           else
           {
               // Count as FULL cycle
               cycle_process_counts(...);
               // Remove A and B (two points)
               residue_cnt -= 2;
           }
       }
       else break;
   }

**Algorithm Steps:**

1. Read the last three turning points from residue: A, B, C
2. Calculate ranges: Y = |A - B|, X = |B - C|
3. Check closure condition: ``X >= Y``
4. If closed:

   - Check if range Y includes the first turning point Z
   - If Y includes Z: count as **half cycle**, remove only A
   - If Y does not include Z: count as **full cycle**, remove A and B
   - Repeat from step 1

5. If not closed, stop checking (new turning point needed)

**Special Handling of First Point (Z):**

The ASTM method has unique handling for cycles that "wrap around" to include
the first turning point of the measurement:

- If the cycle range includes Z, it's counted as a **half cycle** (0.5)
- This prevents double-counting when the signal is repeated or periodic

Key Differences
================

Comparison Table
----------------

+------------------------+---------------------------+---------------------------+
| Aspect                 | 4-Point (DIN 45667)       | 3-Point (ASTM E 1049)     |
+========================+===========================+===========================+
| Points examined        | 4 (A, B, C, D)            | 3 (A, B, C)               |
+------------------------+---------------------------+---------------------------+
| Closure condition      | min(A,D) <= min(B,C) &&   | \|B-C\| >= \|A-B\|        |
|                        | max(B,C) <= max(A,D)      |                           |
+------------------------+---------------------------+---------------------------+
| Cycle counting         | Always full cycles        | Full or half cycles       |
+------------------------+---------------------------+---------------------------+
| First point handling   | No special handling       | Half-cycle if included    |
+------------------------+---------------------------+---------------------------+
| Points removed         | Always 2 (B and C)        | 1 or 2 depending on Z     |
+------------------------+---------------------------+---------------------------+
| Minimum residue        | 4 points                  | 3 points                  |
+------------------------+---------------------------+---------------------------+

Half-Cycle Counting (ASTM Only)
--------------------------------

The ASTM method counts cycles as half-cycles (0.5) when the cycle range includes
the first turning point of the measurement. This is intended to handle edge effects
at measurement boundaries.

**When half-cycles occur:**

A cycle A-B is counted as 0.5 instead of 1.0 when the range [min(A,B), max(A,B)]
contains the value of the first turning point Z.

DIN 45667 with Repeated Residue
--------------------------------

The 4-point method (DIN 45667) does not have built-in half-cycle counting.
However, it achieves correct edge handling through a different mechanism:
the **Repeated residue method** (``RFC_RES_REPEATED``).

When using ``RFC_RES_REPEATED``, the algorithm:

1. Appends the residue to itself (simulating signal repetition)
2. Extracts additional **full cycles** from the doubled residue
3. These cycles are counted as full cycles (weight 1.0)

This approach is **intuitive and physically meaningful**: The results are
exactly equivalent to actually repeating the time series. This makes the
Repeated method ideal for **extrapolation** scenarios where a measured load
sequence represents a repeating pattern (e.g., one revolution, one duty cycle).

**Comparison with ASTM:**

- **ASTM**: Counts affected cycles as 0.5 during normal processing
- **DIN + Repeated**: Counts all extracted cycles as full cycles (1.0)

The key difference: DIN 45667 with ``RFC_RES_REPEATED`` counts cycles from
the residue as full cycles, while ASTM counts edge-affected cycles as half-cycles.
This leads to higher damage estimates with DIN 45667.

**Extrapolation Behavior:**

When extrapolating results (multiplying cycle counts for longer service life),
the methods behave very differently:

- **DIN + Repeated**: Scales correctly. If the time series is repeated N times,
  the damage scales by factor N. The residue handling produces more accurate
  results.

- **ASTM**: Underestimates damage increasingly with higher extrapolation factors.
  The half-cycles at the boundaries remain at 0.5 regardless of repetition count,
  while physically they should contribute more damage with each repetition.

**Example:**

Consider a measurement representing one machine cycle, extrapolated to 1000 cycles:

- **DIN + Repeated**: Residue cycles counted as 1.0 × 1000 = 1000 cycles
- **ASTM**: Edge cycles counted as 0.5, remaining 0.5 even after extrapolation

The ASTM underestimation grows proportionally with the extrapolation factor.

**Damage Estimation Comparison:**

+------------------------+---------------------------+---------------------------+
| Scenario               | ASTM 3-Point              | DIN 4-Point + Repeated    |
+========================+===========================+===========================+
| Edge-affected cycles   | Counted as 0.5            | Full cycles (1.0)         |
+------------------------+---------------------------+---------------------------+
| Residue cycles         | Some counted as 0.5       | All counted as 1.0        |
+------------------------+---------------------------+---------------------------+
| Extrapolation          | Underestimates            | Scales correctly          |
|                        | (error grows with factor) | (physically accurate)     |
+------------------------+---------------------------+---------------------------+
| Total damage           | Tends to underestimate    | Slightly conservative     |
+------------------------+---------------------------+---------------------------+
| Safety margin          | Lower                     | Higher (preferred)        |
+------------------------+---------------------------+---------------------------+

For fatigue analysis, a slightly conservative result (higher damage estimate)
is generally preferred over underestimation, making the DIN 45667 method with
``RFC_RES_REPEATED`` the safer choice for critical applications.

Algorithmic Complexity
-----------------------

Both algorithms have similar complexity:

- **Time Complexity:** O(n) for n turning points

  - Each turning point is processed once
  - Cycle extraction may require multiple iterations, but total work is O(n)

- **Space Complexity:** O(r) where r is residue size

  - Both methods maintain a residue buffer
  - Residue size bounded by 2 × class_count in worst case

The 4-point method is **not** O(n²) as I previously stated incorrectly.
Both methods process the residue from the end without backtracking.

Result Equivalence
==================

Same Cycles, Different Counting
--------------------------------

For signals where no cycle includes the first turning point, both methods
produce **identical results**.

The only difference occurs when a cycle's range includes the first point:

- 4-Point: Counts as 1.0 (full cycle)
- ASTM: Counts as 0.5 (half cycle)

**Example where results differ:**

.. code-block:: text

   Signal: [5, 10, 0, 15, 5]
           Z=5 (first point)

   Cycle found: (10, 0) with range 10

   Range [0, 10] includes Z=5

   4-Point result: 1.0 full cycle
   ASTM result:    0.5 half cycle

For most practical signals (especially long measurements), the difference
is negligible since only cycles near the measurement boundaries are affected.

**Important Note on Damage Estimation:**

When comparing total damage results:

- **ASTM tends to underestimate** damage due to half-cycle counting
- **DIN 45667 with RFC_RES_REPEATED** provides a more conservative estimate

For safety-critical applications, the slightly higher damage estimate from
DIN 45667 + Repeated residue is preferred, as it provides a safety margin
against fatigue failure.

Verification
------------

.. code-block:: c

   rfc_ctx_s ctx_4pt, ctx_astm;

   // Initialize both with same parameters
   RFC_init(&ctx_4pt, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
   RFC_init(&ctx_astm, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Set counting methods
   ctx_4pt.counting_method = RFC_COUNTING_METHOD_4PTM;
   ctx_astm.counting_method = RFC_COUNTING_METHOD_ASTM;

   // Feed same data
   RFC_feed(&ctx_4pt, data, data_len);
   RFC_feed(&ctx_astm, data, data_len);

   // Finalize
   RFC_finalize(&ctx_4pt, RFC_RES_IGNORE);
   RFC_finalize(&ctx_astm, RFC_RES_IGNORE);

   // Compare - may differ slightly due to half-cycle counting
   printf("4-Point: full=%.1f, half=%.1f\n",
          ctx_4pt.full_inc, ctx_4pt.half_inc);
   printf("ASTM:    full=%.1f, half=%.1f\n",
          ctx_astm.full_inc, ctx_astm.half_inc);

When to Use Each Method
=======================

Use 4-Point Method When:
-------------------------

- **Simplicity is preferred** - No half-cycle complexity
- **DIN 45667 compliance** - German standard requirement
- **FVA-Richtlinie compliance** - German drive technology guideline
- **All cycles should be full** - No half-cycle counting desired
- **Conservative damage estimate** - Use with ``RFC_RES_REPEATED`` for safety
- **Default choice** - This is the library default

Use ASTM Method When:
---------------------

- **ASTM E 1049 compliance** - US standard requirement
- **Half-cycle distinction needed** - Separate full/half cycle counts

**Recommendation:**

For safety-critical fatigue analysis, use the **4-point method (DIN 45667)**
with ``RFC_RES_REPEATED`` residue handling. This combination provides:

- Correct handling of edge/boundary cycles
- Slightly conservative damage estimates
- No risk of underestimating fatigue damage

Implementation in This Library
===============================

Enabling Methods
----------------

**Compile-Time (ASTM support is optional):**

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_ASTM_SUPPORT=ON

**C Code:**

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Use 4-point method (default)
   ctx.counting_method = RFC_COUNTING_METHOD_4PTM;

   // Or use ASTM method
   ctx.counting_method = RFC_COUNTING_METHOD_ASTM;

   RFC_feed(&ctx, data, data_len);
   RFC_finalize(&ctx, RFC_RES_REPEATED);

**C++ Code:**

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0, 0.0, 1.0);

   // Use ASTM method
   rf.ctx_get().counting_method = RF::RFC_COUNTING_METHOD_ASTM;

   rf.feed(data, data_len);
   rf.finalize(Rainflow::RFC_RES_REPEATED);

Available Methods
-----------------

.. code-block:: c

   typedef enum rfc_counting_method
   {
       RFC_COUNTING_METHOD_NONE       = 0,  // No counting
       RFC_COUNTING_METHOD_4PTM       = 1,  // 4-point method (default)
       RFC_COUNTING_METHOD_HCM        = 2,  // HCM method (Clormann/Seeger)
       RFC_COUNTING_METHOD_ASTM       = 3,  // ASTM 3-point method
       RFC_COUNTING_METHOD_DELEGATED  = 4,  // Custom delegate function
   } rfc_counting_method_e;

Standards Compliance
====================

4-Point Method Complies With:
------------------------------

- **DIN 45667** - German standard for load spectrum classification
- **FVA-Richtlinie** - German drive technology association guideline
- Most European fatigue analysis standards

ASTM Method Complies With:
---------------------------

- **ASTM E 1049-85 (2011)** - US standard for cycle counting
- SAE fatigue analysis recommendations
- Many aerospace industry requirements

See Also
========

- `algorithm.rst <algorithm.rst>`_ - Overview of rainflow counting
- `features.rst <features.rst>`_ - Compile-time feature selection
- `residue_methods.rst <residue_methods.rst>`_ - Handling unclosed cycles
- `references.rst <references.rst>`_ - ASTM E 1049 and DIN 45667 citations

References
==========

**ASTM E 1049-85 (2011)**
   "Standard Practices for Cycle Counting in Fatigue Analysis"
   ASTM International, West Conshohocken, PA

**DIN 45667:2010**
   "Classification counting methods for one-parameter and two-parameter random loadings"

**FVA-Richtlinie (2010)**
   "Zählverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"
   Forschungsvereinigung Antriebstechnik e.V.

**Clormann & Seeger (1985)**
   "Rainflow-HCM / Ein Hysteresisschleifen-Zählalgorithmus auf werkstoffmechanischer Grundlage"
   TU Darmstadt

For implementation details:

- ``rainflow.c`` (``cycle_find_4ptm()``)
- ``rainflow.c`` (``cycle_find_astm()``)
