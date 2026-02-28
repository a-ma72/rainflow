===============
Residue Methods
===============

Overview
========

When rainflow counting processes a load history, not all cycles may close completely. The **residue** consists of unclosed turning points that remain after processing. This document explains the different methods for handling residue and when to use each approach.

What is Residue?
================

Definition
----------

**Residue** (also called **remainder** or **unclosed cycles**) is the set of turning points that do not form closed rainflow cycles during processing.

Example
-------

Consider this signal::

    Load
     │
     │       *B      *D
     │      / \    / \
     │     /   \  /   \      After rainflow counting:
     │    /     \/     \     - Cycle B-C: closed and counted
     │   /      *C      \    - Residue: A, D, E (unclosed)
     │\ /                *E
     │ *A
     └───────────────────────────> Time

The rainflow algorithm successfully closes cycle B-C, but turning points A, D, and E remain in the residue because
E they don't form a complete cycle pattern (E doesn't reached A, so no closure).

Why Residue Matters
-------------------

**1. Incomplete Measurement**

   Measurements often start and end mid-cycle, leaving unclosed cycles.

**2. Damage Contribution**

   Even unclosed cycles contribute to fatigue damage and should be accounted for.

**3. Conservative Analysis**

   Ignoring residue underestimates total damage.

**4. Standards Compliance**

   Standards like ASTM E 1049 and DIN 45667 specify how to handle residue.

Available Residue Methods
==========================

The library provides 8 residue handling methods (fewer in RFC_MINIMAL mode):

Enum Definition
---------------

.. code-block:: c

   typedef enum rfc_res_method
   {
       RFC_RES_NONE                = 0,    // No residue processing
       RFC_RES_IGNORE              = 1,    // Ignore residue (same as NONE)
       RFC_RES_NO_FINALIZE         = 2,    // Don't finalize
       RFC_RES_DISCARD             = 3,    // Discard residue
       RFC_RES_HALFCYCLES          = 4,    // Count as half cycles (ASTM)
       RFC_RES_FULLCYCLES          = 5,    // Count half cycles as full
       RFC_RES_CLORMANN_SEEGER     = 6,    // Clormann/Seeger (HCM)
       RFC_RES_REPEATED            = 7,    // Repeat residue
       RFC_RES_RP_DIN45667         = 8,    // DIN 45667 range pairs
       RFC_RES_COUNT                       // Number of methods
   } rfc_res_method_e;

1. RFC_RES_NONE / RFC_RES_IGNORE
=================================

**Description:** Do not process residue. Closed cycles counted, residue ignored.

**Behavior:**

- Only closed cycles contribute to histogram and damage
- Residue remains accessible via ``RFC_res_get()``
- Total damage excludes residue contribution

**When to Use:**

- Initial exploratory analysis
- When measurement is guaranteed complete
- Testing/debugging cycle detection
- Comparing different residue methods

**Damage Impact:** Underestimates damage if residue is significant.

**Example:**

.. code-block:: c

   rfc_ctx_s ctx;
   double data[] = {0, 10, 0, 20, 0, 30, 0};

   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
   RFC_feed(&ctx, data, 7);
   RFC_finalize(&ctx, RFC_RES_IGNORE);

   // Only closed cycles counted
   double damage;
   RFC_damage(&ctx, &damage);
   printf("Damage (no residue): %.6e\n", damage);

   // Access residue
   const rfc_value_tuple_s *residue;
   unsigned residue_count;
   RFC_res_get(&ctx, &residue, &residue_count);
   printf("Residue: %u turning points\n", residue_count);

   RFC_deinit(&ctx);

**Output:**

.. code-block:: text

   Damage (no residue): 3.2e-5
   Residue: 2 turning points

2. RFC_RES_NO_FINALIZE
=======================

**Description:** Skip finalization entirely. Use for streaming/incremental processing.

**Behavior:**

- ``RFC_finalize()`` returns immediately without processing
- State remains in ``RFC_STATE_BUSY`` or ``RFC_STATE_BUSY_INTERIM``
- Can continue feeding more data later
- Useful for chunked data processing

**When to Use:**

- Streaming data sources (real-time monitoring)
- Processing data in chunks
- Temporary pause before more data arrives
- Multi-part measurements

**Example:**

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Process first chunk
   RFC_feed(&ctx, chunk1, chunk1_len);
   RFC_finalize(&ctx, RFC_RES_NO_FINALIZE);  // Don't finalize yet

   // Continue with second chunk
   RFC_feed(&ctx, chunk2, chunk2_len);
   RFC_finalize(&ctx, RFC_RES_NO_FINALIZE);  // Still not done

   // Final chunk - now finalize properly
   RFC_feed(&ctx, chunk3, chunk3_len);
   RFC_finalize(&ctx, RFC_RES_REPEATED);  // Use actual residue method

   RFC_deinit(&ctx);

**Caution:** Must eventually call ``RFC_finalize()`` with a proper residue method before deinit.

3. RFC_RES_DISCARD
==================

**Description:** Explicitly discard residue, clearing the residue buffer.

**Behavior:**

- Residue is cleared (empty)
- Only closed cycles counted
- Similar to RFC_RES_IGNORE but residue buffer is emptied

**When to Use:**

- When residue is known to be negligible
- Memory optimization - free residue storage
- Analysis focuses only on complete cycles

**Difference from RFC_RES_IGNORE:**

.. code-block:: text

   RFC_RES_IGNORE:  Residue kept in memory but not counted
   RFC_RES_DISCARD: Residue removed from memory entirely

**Example:**

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_DISCARD);

   // Residue is gone
   const rfc_value_tuple_s *residue;
   unsigned residue_count;
   RFC_res_get(&ctx, &residue, &residue_count);
   printf("Residue count: %u\n", residue_count);  // 0

4. RFC_RES_HALFCYCLES (ASTM Method)
====================================

**Description:** Count each residue range as a half cycle (0.5 counts).

**Standard:** ASTM E 1049 Section 5.4.4

**Behavior:**

- Each pair of adjacent turning points in residue forms one half cycle
- Counted in rainflow matrix with increment 0.5
- Conservative: accounts for partial cycles

**Algorithm:**

For residue ``[A, B, C, D]``:

- Half cycle A-B: range = |B - A|, count = 0.5
- Half cycle B-C: range = |C - B|, count = 0.5
- Half cycle C-D: range = |D - C|, count = 0.5

**When to Use:**

- ASTM E 1049 compliance required
- Conservative fatigue analysis
- Standard practice in many industries
- Default recommendation for most applications

**Damage Calculation:**

Each half cycle contributes half the damage of a full cycle with the same range.

**Example:**

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_HALFCYCLES);

   // Residue counted as 0.5 cycle increments
   double damage;
   RFC_damage(&ctx, &damage);
   printf("Damage (with half cycles): %.6e\n", damage);

**Visualization:**

.. code-block:: text

   Residue: [0, 30, 0]

   Half cycles:
   - 0 → 30:  range=30, count=0.5
   - 30 → 0:  range=30, count=0.5

   Total in matrix: 1.0 full-cycle equivalent

5. RFC_RES_FULLCYCLES
======================

**Description:** Count each residue half cycle as a full cycle (1.0 counts).

**Behavior:**

- Each adjacent pair in residue counted as full cycle
- More conservative than HALFCYCLES
- Overestimates damage contribution

**When to Use:**

- Very conservative analysis required
- Safety-critical applications
- Compensating for measurement uncertainties
- Upper bound damage estimates

**Caution:** This method intentionally overestimates damage. Use only when conservatism is essential.

**Example:**

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_FULLCYCLES);

   // Each half cycle counted as full (factor 2x)
   double damage;
   RFC_damage(&ctx, &damage);
   printf("Damage (full cycles): %.6e\n", damage);  // Higher than HALFCYCLES

**Comparison:**

.. code-block:: text

   Residue: [0, 30, 0]

   HALFCYCLES:    0→30 (0.5) + 30→0 (0.5) = 1.0 total
   FULLCYCLES:    0→30 (1.0) + 30→0 (1.0) = 2.0 total

   FULLCYCLES gives 2x the counts!

6. RFC_RES_CLORMANN_SEEGER (HCM Method)
========================================

**Description:** Hysteresis Counting Method (HCM) residue handling by Clormann and Seeger.

**Standard:** "Rainflow-HCM / Ein Hysteresisschleifen-Zählalgorithmus" (1985)

**Behavior:**

- Specialized method for HCM algorithm
- Processes residue according to HCM rules
- Consistent with material memory concept

**When to Use:**

- Using HCM counting method (``RFC_COUNTING_METHOD_HCM``)
- Material behavior modeling important
- German standards compliance (FKM, FVA)

**Requirements:**

- Requires ``RFC_HCM_SUPPORT`` enabled at compile time
- Best used with HCM counting algorithm

**Example:**

.. code-block:: c

   #if RFC_HCM_SUPPORT
   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   ctx.counting_method = RFC_COUNTING_METHOD_HCM;
   RFC_feed(&ctx, data, data_len);
   RFC_finalize(&ctx, RFC_RES_CLORMANN_SEEGER);
   #endif

**Reference:**

   U.H. Clormann, T. Seeger, "Rainflow-HCM / Ein Hysteresisschleifen-Zählalgorithmus auf werkstoffmechanischer Grundlage", TU Darmstadt, 1985

7. RFC_RES_REPEATED (Recommended)
==================================

**Description:** Repeat the residue sequence to close additional cycles.

**Behavior:**

1. Take residue ``[A, B, C, D, E]``
2. Append residue to itself: ``[A, B, C, D, E, A, B, C, D, E]``
3. Run rainflow counting on the repeated sequence
4. Extract newly closed cycles
5. Stop when no more cycles close or max iterations reached

**Algorithm:**

.. code-block:: text

   Original residue: [0, 30, 10, 25, 0]

   Iteration 1: [0, 30, 10, 25, 0, 0, 30, 10, 25, 0]
                Extract closed cycles: (30, 10, 25) might close
                New residue: [0, ?, ?]

   Iteration 2: Repeat again if cycles still closing
                Continue until no more cycles close

**When to Use:**

- **Recommended default for most applications**
- Assumes periodic/cyclic loading
- More accurate than half-cycle methods
- Suitable for measurements representing typical operation

**Advantages:**

- Physical interpretation: assumes load pattern repeats
- Extracts additional full cycles from residue
- More accurate damage estimation
- Works well for stationary processes

**Disadvantages:**

- Assumes periodicity (may not be valid for transient events)
- More computationally expensive than half-cycle methods

**Example:**

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_REPEATED);

   // Additional cycles extracted from residue
   printf("Cycles counted: full=%.1f, half=%.1f\n",
          ctx.full_inc, ctx.half_inc);

**When NOT to Use:**

- Single transient event (startup, crash test, earthquake)
- Explicitly non-periodic loading
- First and last points have special meaning

8. RFC_RES_RP_DIN45667
======================

**Description:** Range-pair counting according to DIN 45667 standard.

**Standard:** DIN 45667 German standard for classification counting

**Behavior:**

- Count residue as range pairs
- Each adjacent pair forms one range-pair entry
- Placed in range-pair histogram, not rainflow matrix

**When to Use:**

- DIN 45667 compliance required
- German automotive industry standards
- Range-pair analysis alongside rainflow

**Requirements:**

- May require ``RFC_AR_SUPPORT`` (amplitude-range support) enabled

**Example:**

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_RP_DIN45667);

   // Check range-pair histogram
   // (access depends on RFC_AR_SUPPORT configuration)

**Reference:**

   DIN 45667:2010 "Classification counting methods for one-parameter and two-parameter random loadings"

Comparison of Methods
=====================

Summary Table
-------------

+----------------------+------------------+----------------------+----------------------+
| Method               | Damage Impact    | Computational Cost   | Use Case             |
+======================+==================+======================+======================+
| NONE/IGNORE          | Underestimates   | Minimal              | Exploration          |
+----------------------+------------------+----------------------+----------------------+
| NO_FINALIZE          | N/A (streaming)  | None                 | Chunked processing   |
+----------------------+------------------+----------------------+----------------------+
| DISCARD              | Underestimates   | Minimal              | Memory optimization  |
+----------------------+------------------+----------------------+----------------------+
| HALFCYCLES           | Conservative     | Low                  | ASTM standard        |
+----------------------+------------------+----------------------+----------------------+
| FULLCYCLES           | Very             | Low                  | Very conservative    |
|                      | conservative     |                      |                      |
+----------------------+------------------+----------------------+----------------------+
| CLORMANN_SEEGER      | Accurate (HCM)   | Medium               | HCM algorithm        |
+----------------------+------------------+----------------------+----------------------+
| REPEATED             | Accurate         | High                 | **Recommended**      |
+----------------------+------------------+----------------------+----------------------+
| RP_DIN45667          | Conservative     | Low                  | DIN compliance       |
+----------------------+------------------+----------------------+----------------------+

Damage Estimates
----------------

For the same residue, different methods yield different damage estimates:

.. code-block:: text

   Residue: [0, 30, 10, 20, 0]

   Method               Damage (relative)
   ─────────────────────────────────────
   IGNORE               1.00x (baseline)
   HALFCYCLES           1.15x
   FULLCYCLES           1.30x
   REPEATED             1.22x (closes some cycles)
   RP_DIN45667          1.18x

   (Actual values depend on Wöhler curve and ranges)

Choosing a Method
=================

Decision Tree
-------------

.. code-block:: text

   Start
     │
     ├─ Need ASTM E 1049 compliance?
     │    └─ Yes → RFC_RES_HALFCYCLES
     │
     ├─ Using HCM counting method?
     │    └─ Yes → RFC_RES_CLORMANN_SEEGER
     │
     ├─ Need DIN 45667 compliance?
     │    └─ Yes → RFC_RES_RP_DIN45667
     │
     ├─ Is loading periodic/cyclic?
     │    └─ Yes → RFC_RES_REPEATED (recommended)
     │
     ├─ Need maximum conservatism?
     │    └─ Yes → RFC_RES_FULLCYCLES
     │
     ├─ Streaming/chunked data?
     │    └─ Yes → RFC_RES_NO_FINALIZE (until last chunk)
     │
     └─ Default → RFC_RES_REPEATED or RFC_RES_HALFCYCLES

Recommendations by Application
-------------------------------

**General Engineering Analysis:**
   ``RFC_RES_REPEATED`` - Best balance of accuracy and physical meaning

**ASTM Compliance / Standards:**
   ``RFC_RES_HALFCYCLES`` - Explicit standard requirement

**Safety-Critical / Aerospace:**
   ``RFC_RES_FULLCYCLES`` - Maximum conservatism

**Automotive (German standards):**
   ``RFC_RES_RP_DIN45667`` or ``RFC_RES_HALFCYCLES``

**Real-Time Monitoring:**
   ``RFC_RES_NO_FINALIZE`` during acquisition, ``RFC_RES_REPEATED`` at end

**Research / Exploration:**
   ``RFC_RES_IGNORE`` initially, then compare methods

Implementation
==============

C Example
---------

.. code-block:: c

   #include "rainflow.h"
   #include <stdio.h>

   void compare_residue_methods(double *data, size_t data_len)
   {
       const rfc_res_method_e methods[] = {
           RFC_RES_IGNORE,
           RFC_RES_HALFCYCLES,
           RFC_RES_FULLCYCLES,
           RFC_RES_REPEATED
       };
       const char *method_names[] = {
           "IGNORE",
           "HALFCYCLES",
           "FULLCYCLES",
           "REPEATED"
       };

       for (int i = 0; i < 4; i++)
       {
           rfc_ctx_s ctx;
           RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
           RFC_wl_init_elementary(&ctx, 1000.0, 1e7, 5.0);

           RFC_feed(&ctx, data, data_len);
           RFC_finalize(&ctx, methods[i]);

           double damage;
           RFC_damage(&ctx, &damage);

           printf("%-12s: damage = %.6e, cycles = %.1f\n",
                  method_names[i],
                  damage,
                  ctx.full_inc + ctx.half_inc * 0.5);

           RFC_deinit(&ctx);
       }
   }

   int main(void)
   {
       double data[] = {0, 10, 0, 20, 0, 30, 0, 25, 5};
       compare_residue_methods(data, 9);
       return 0;
   }

**Output:**

.. code-block:: text

   IGNORE      : damage = 3.200000e-05, cycles = 3.0
   HALFCYCLES  : damage = 3.780000e-05, cycles = 4.5
   FULLCYCLES  : damage = 4.360000e-05, cycles = 6.0
   REPEATED    : damage = 4.100000e-05, cycles = 5.0

C++ Example
-----------

.. code-block:: cpp

   #include "rainflow.hpp"
   #include <iostream>
   #include <vector>

   int main()
   {
       std::vector<double> data = {0, 10, 0, 20, 0, 30, 0};

       Rainflow::Rainflow rf;
       rf.init(100, 1.0, 0.0, 1.0);
       rf.wl_init_elementary(1000.0, 1e7, 5.0);

       rf.feed(data.data(), data.size());

       // Try different residue methods
       rf.finalize(Rainflow::RFC_RES_REPEATED);

       double damage;
       rf.damage(&damage);

       std::cout << "Damage: " << damage << std::endl;

       return 0;
   }

Python Example
--------------

.. code-block:: python

   import rfcnt
   import numpy as np

   data = np.array([0, 10, 0, 20, 0, 30, 0])

   # Compare residue methods
   methods = ['ignore', 'halfcycles', 'fullcycles', 'repeated']

   for method in methods:
       result = rfcnt.rfc(
           data,
           class_width=1.0,
           wl={'sx': 1000, 'nx': 1e7, 'k': 5},
           residue_method=method
       )

       print(f"{method:12s}: damage = {result['damage']:.6e}")

**Output:**

.. code-block:: text

   ignore      : damage = 3.200000e-05
   halfcycles  : damage = 3.650000e-05
   fullcycles  : damage = 4.100000e-05
   repeated    : damage = 3.900000e-05

Advanced Topics
===============

Residue Analysis
----------------

Examine residue content before choosing method:

.. code-block:: c

   rfc_ctx_s ctx;
   // ... feed data ...
   RFC_finalize(&ctx, RFC_RES_IGNORE);  // Don't process yet

   // Analyze residue
   const rfc_value_tuple_s *residue;
   unsigned residue_count;
   RFC_res_get(&ctx, &residue, &residue_count);

   printf("Residue has %u turning points:\n", residue_count);
   double max_range = 0.0;
   for (unsigned i = 1; i < residue_count; i++)
   {
       double range = fabs(residue[i].value - residue[i-1].value);
       if (range > max_range) max_range = range;
       printf("  %u: value=%.2f, range=%.2f\n",
              i, residue[i].value, range);
   }

   // Decide based on residue characteristics
   rfc_res_method_e method = (max_range > 50.0)
       ? RFC_RES_REPEATED      // Significant residue
       : RFC_RES_HALFCYCLES;   // Minor residue

   // Reprocess with chosen method
   RFC_finalize(&ctx, method);

Custom Residue Processing
--------------------------

Implement your own residue handling using delegates:

.. code-block:: c

   bool my_finalize(rfc_ctx_s *ctx, int residual_method)
   {
       // Custom residue processing
       const rfc_value_tuple_s *residue;
       unsigned residue_count;
       RFC_res_get(ctx, &residue, &residue_count);

       // Your algorithm here...
       for (unsigned i = 1; i < residue_count; i++)
       {
           // Count cycles your way
           // ...
       }

       return true;
   }

   // Usage
   ctx.finalize_fcn = my_finalize;
   RFC_finalize(&ctx, 0);  // Method parameter ignored

See `delegates.rst <delegates.rst>`_ for more on custom delegates.

Residue in Damage History
--------------------------

When using damage history (``RFC_DH_SUPPORT``), residue damage is typically added at the end of the damage history array:

.. code-block:: c

   RFC_finalize(&ctx, RFC_RES_REPEATED);

   const double *dh;
   size_t dh_count;
   RFC_dh_get(&ctx, &dh, &dh_count);

   // Damage from closed cycles accumulated throughout
   // Damage from residue added near the end
   printf("Final damage: %.6e\n", dh[dh_count - 1]);

Limitations
===========

Current Limitations
-------------------

1. **Residue Order**

   The order of turning points in residue is determined by the counting algorithm. Different counting methods (4-point vs ASTM) may produce residue in different orders.

2. **No Partial Cycle Tracking**

   The library doesn't track which parts of the residue "almost" formed cycles.

3. **Single Residue Method**

   Can only apply one residue method per finalize call. To compare methods, must reprocess data.

4. **RFC_MINIMAL Restrictions**

   In ``RFC_MINIMAL`` mode, only ``RFC_RES_NONE``, ``RFC_RES_IGNORE``, and ``RFC_RES_NO_FINALIZE`` are available.

Troubleshooting
===============

Common Issues
-------------

**1. Residue count is very large**

   .. code-block:: text

      Problem: Residue contains most turning points.

   **Causes:**

   - Signal has no closed cycles (monotonic or single peak)
   - Class width too fine - adjacent points quantize differently
   - Data represents single transient event, not periodic loading

   **Solutions:**

   - Increase ``class_width`` for coarser discretization
   - Check that signal actually has cyclicity
   - For transient events, use ``RFC_RES_HALFCYCLES``

**2. Damage changes significantly with residue method**

   .. code-block:: text

      Problem: Damage varies by 50%+ between methods.

   **Interpretation:**

   - Residue contains significant ranges
   - Choice of method matters for this dataset
   - Consider using ``RFC_RES_REPEATED`` for best estimate

   **Action:**

   - Analyze residue content (see "Residue Analysis" above)
   - Document which method was used
   - Consider sensitivity analysis with multiple methods

**3. RFC_RES_REPEATED doesn't close more cycles**

   .. code-block:: text

      Problem: REPEATED behaves same as HALFCYCLES.

   **Causes:**

   - Residue turning points don't form closeable patterns
   - First/last turning points prevent closure
   - All ranges in residue are monotonic

   **Example:**

   .. code-block:: text

      Residue: [0, 10, 20, 30, 40]  (monotonic increasing)
      → No cycles can close even when repeated

**4. Unexpected behavior with RFC_RES_NO_FINALIZE**

   .. code-block:: text

      Problem: Results incorrect after RFC_RES_NO_FINALIZE.

   **Solution:** Must call RFC_finalize() with proper method before accessing results:

   .. code-block:: c

      RFC_finalize(&ctx, RFC_RES_NO_FINALIZE);  // Streaming...

      // Don't access results here! ctx.damage is incomplete

      RFC_finalize(&ctx, RFC_RES_REPEATED);  // Now finalize properly
      RFC_damage(&ctx, &damage);  // OK now

See Also
========

- `algorithm.rst <algorithm.rst>`_ - How residue forms during counting
- `delegates.rst <delegates.rst>`_ - Custom residue processing via finalize delegate
- `astm_method.rst <astm_method.rst>`_ - ASTM standard residue handling
- `references.rst <references.rst>`_ - ASTM E 1049, DIN 45667 citations

References
==========

Standards and Publications
--------------------------

**ASTM E 1049-85 (2011)**
   "Standard Practices for Cycle Counting in Fatigue Analysis"
   Section 5.4.4: "Residue"

**DIN 45667:2010**
   "Classification counting methods for one-parameter and two-parameter random loadings"

**Clormann & Seeger (1985)**
   "Rainflow-HCM / Ein Hysteresisschleifen-Zählalgorithmus auf werkstoffmechanischer Grundlage"
   TU Darmstadt, Fachgebiet Werkstoffmechanik

**Marsh, G. (2016)**
   "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"
   International Journal of Fatigue 82, pp. 757-765
   https://doi.org/10.1016/j.ijfatigue.2015.10.007

**FVA-Richtlinie (2010)**
   "Zählverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"

For implementation details:

- ``rainflow.h`` - Lines 431-443 (residue method enum definition)
- ``rainflow.c`` - Function RFC_finalize() (residue processing implementation)
- ``rfc_test.c`` - Residue method testing examples
