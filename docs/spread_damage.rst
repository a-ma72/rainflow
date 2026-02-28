==============
Spread Damage
==============

Overview
========

When tracking damage history, each closed rainflow cycle contributes damage at specific time points. **Spread damage** refers to how this cycle damage is distributed across the turning points involved in the cycle. Different spreading methods represent different physical interpretations of when and how damage accumulates.

This feature is only available with ``RFC_DH_SUPPORT`` enabled and is used in conjunction with damage history tracking.

What is Spread Damage?
=======================

The Problem
-----------

When a rainflow cycle closes (e.g., cycle B-C in the 4-point algorithm), it has a computed damage value D. But **when** did this damage occur?

.. code-block:: text

                     * D (t=4)
                    / \
             B *<--/          Cycle B-C closes when D arrives
              / \ /           Damage for B-C = 0.001
             /   * C (t=3)
          \ /
           * A (t=1)    (t=2)

**Question:** Should the damage 0.001 be assigned to:

- Time t=2 (when B occurred)?
- Time t=3 (when C occurred)?
- Spread between t=2 and t=3?
- Spread between t=2 and t=4 (when cycle closed)?

Different spreading methods answer this question differently.

Physical Interpretation
-----------------------

The choice of spreading method reflects assumptions about the fatigue process:

**Instantaneous Models** (FULL_P2, FULL_P3)
   Damage occurs at a single instant (peak or valley)

**Linear Models** (RAMP variants)
   Damage accumulates linearly over time or amplitude

**Transient Models** (TRANSIENT_23, TRANSIENT_23c)
   Damage follows a transient response based on material behavior

Available Spread Damage Methods
================================

The library provides 10 spreading methods (9 methods plus NONE):

Enum Definition
---------------

.. code-block:: c

   typedef enum rfc_sd_method
   {
       RFC_SD_NONE                 = -1,   // No damage spreading
       RFC_SD_HALF_23              =  0,   // Split 50/50 between P2 and P3
       RFC_SD_RAMP_AMPLITUDE_23    =  1,   // Ramp by amplitude P2→P3
       RFC_SD_RAMP_DAMAGE_23       =  2,   // Linear ramp P2→P3
       RFC_SD_RAMP_AMPLITUDE_24    =  3,   // Ramp by amplitude P2→P4
       RFC_SD_RAMP_DAMAGE_24       =  4,   // Linear ramp P2→P4
       RFC_SD_FULL_P2              =  5,   // All damage at P2
       RFC_SD_FULL_P3              =  6,   // All damage at P3
       RFC_SD_TRANSIENT_23         =  7,   // Transient P2→P3
       RFC_SD_TRANSIENT_23c        =  8,   // Transient P2→P3c (recommended)
       RFC_SD_COUNT                         // Number of methods
   } rfc_sd_method_e;

Rainflow Cycle Points Reference
--------------------------------

In the 4-point rainflow algorithm, a cycle is defined by 4 consecutive turning points:

.. code-block:: text

                     * P4 (D in notation)
                    / \
            P2 *<--/          Cycle P2-P3 closes when P4 arrives
              / \ /           P1-P2-P3-P4 are examined
             /   * P3 (C)
          \ /
           * P1 (A)

**Points:**

- **P1**: First point (before cycle)
- **P2**: Cycle start (peak or valley)
- **P3**: Cycle end (valley or peak)
- **P4**: Point that causes cycle P2-P3 to close
- **P3c**: The position where the cycle actually closes (between P3 and P4)

Method Descriptions
===================

1. RFC_SD_NONE
--------------

**Description:** No damage spreading. Damage history not updated.

**Behavior:**

- Damage is calculated and accumulated in total
- Damage history array remains zero or unchanged
- Minimal computational cost

**When to Use:**

- Only total damage needed, not damage history
- Performance optimization when history not required

**Example:**

.. code-block:: c

   RFC_dh_init(&ctx, RFC_SD_NONE, NULL, 0, false);
   // Damage history disabled

2. RFC_SD_HALF_23
-----------------

**Description:** Split cycle damage equally (50/50) between P2 and P3.

**Algorithm:**

.. code-block:: text

   D_total = damage from cycle P2-P3
   D_P2 += D_total / 2
   D_P3 += D_total / 2

**Visualization:**

.. code-block:: text

   Damage
     │
     │    ┌─────────┐
     │    │  D/2    │  D/2
     │────┼─────────┼─────────>
         P2        P3       Time

**Physical Meaning:**

- Damage accumulates half at loading, half at unloading
- Simple, symmetric model
- No consideration of intermediate points

**When to Use:**

- Simple damage tracking
- When physical timing doesn't matter much
- Coarse damage history sufficient

3. RFC_SD_RAMP_AMPLITUDE_23
----------------------------

**Description:** Spread damage proportional to amplitude change between P2 and P3.

**Algorithm:**

For each turning point ``Pi`` between P2 and P3:

.. code-block:: text

   amplitude_ratio = |Pi - P2| / |P3 - P2|
   D_Pi += D_total × amplitude_ratio

**Visualization:**

.. code-block:: text

   Damage
     │       ╱
     │      ╱
     │     ╱
     │    ╱
     │   ╱
     │──┴──────────────>
       P2      P3    Time

   Damage increases linearly with amplitude

**Physical Meaning:**

- Damage proportional to strain/stress magnitude
- More damage near peak amplitude
- Physically motivated for amplitude-dependent damage

**When to Use:**

- Damage is amplitude-sensitive
- Material models with strain dependence
- More physical than RAMP_DAMAGE_23

4. RFC_SD_RAMP_DAMAGE_23
-------------------------

**Description:** Spread damage evenly (linearly in time) between P2 and P3.

**Algorithm:**

For each turning point ``Pi`` between P2 and P3:

.. code-block:: text

   time_ratio = (i - P2_index) / (P3_index - P2_index)
   D_Pi += D_total × time_ratio

**Visualization:**

.. code-block:: text

   Damage
     │       ╱
     │      ╱
     │     ╱
     │    ╱
     │   ╱
     │──┴──────────────>
       P2      P3    Time

   Linear ramp in time (not amplitude)

**Physical Meaning:**

- Damage accumulates uniformly over time
- No preference for high/low amplitude points
- Simplest time-based model

**When to Use:**

- Time-based damage accumulation
- Damage rate is constant
- Simple linear interpolation desired

5. RFC_SD_RAMP_AMPLITUDE_24
----------------------------

**Description:** Spread damage exponentially according to amplitude from P2 to P4.

**Algorithm:**

Spread damage over the full 4-point sequence P2→P3→P4:

.. code-block:: text

   For each point Pi between P2 and P4:
       amplitude_factor = f(|Pi - P2|)  // Exponential or power law
       D_Pi += D_total × amplitude_factor

**Visualization:**

.. code-block:: text

   Damage
     │         ╱╲
     │        ╱  ╲
     │       ╱    ╲
     │      ╱      ╲
     │─────┴────────┴──>
        P2   P3    P4  Time

**Physical Meaning:**

- Accounts for full cycle context (P2 through P4)
- More sophisticated amplitude weighting
- Considers cycle closure dynamics

**When to Use:**

- Detailed damage modeling
- Full cycle context matters
- Complex amplitude-dependent behavior

6. RFC_SD_RAMP_DAMAGE_24
-------------------------

**Description:** Spread damage evenly (linearly) from P2 to P4.

**Algorithm:**

.. code-block:: text

   For each point Pi between P2 and P4:
       time_ratio = (i - P2_index) / (P4_index - P2_index)
       D_Pi += D_total × time_ratio

**Visualization:**

.. code-block:: text

   Damage
     │           ╱
     │          ╱
     │         ╱
     │        ╱
     │───────┴────────>
        P2  P3  P4   Time

**Physical Meaning:**

- Linear accumulation over full cycle duration
- Extends damage to cycle closure point
- Smoother distribution than RAMP_DAMAGE_23

**When to Use:**

- When cycle closure timing matters
- Smoother damage history desired
- Full 4-point context important

7. RFC_SD_FULL_P2
-----------------

**Description:** Assign all cycle damage to P2 (cycle start).

**Algorithm:**

.. code-block:: text

   D_P2 += D_total
   // All other points: no damage

**Visualization:**

.. code-block:: text

   Damage
     │    │
     │    │ D_total
     │────┼──────────>
         P2       Time

   All damage at P2 (impulse)

**Physical Meaning:**

- Damage occurs at maximum stress/strain point
- Loading event causes damage
- Ignores unloading contribution

**When to Use:**

- Damage dominated by peak loading
- Crack initiation models
- Maximum stress criteria

8. RFC_SD_FULL_P3
-----------------

**Description:** Assign all cycle damage to P3 (cycle end).

**Algorithm:**

.. code-block:: text

   D_P3 += D_total
   // All other points: no damage

**Visualization:**

.. code-block:: text

   Damage
     │              │
     │              │ D_total
     │──────────────┼──>
                   P3  Time

**Physical Meaning:**

- Damage occurs at unloading/return point
- Cycle completion triggers damage
- Alternative to FULL_P2

**When to Use:**

- Damage at unloading important
- Cycle closure is critical event
- Symmetric alternative to FULL_P2

9. RFC_SD_TRANSIENT_23 (Advanced)
----------------------------------

**Description:** Spread damage according to transient response between P2 and P3.

**Algorithm:**

Uses an exponential or transient function to model how damage evolves:

.. code-block:: text

   For each point Pi between P2 and P3:
       t_norm = (i - P2_index) / (P3_index - P2_index)
       transient_factor = 1 - exp(-lambda × t_norm)
       D_Pi += D_total × transient_factor

**Visualization:**

.. code-block:: text

   Damage
     │       ╭─────
     │      ╱
     │     ╱
     │    ╱
     │───┴──────────>
       P2      P3   Time

   Exponential buildup (transient response)

**Physical Meaning:**

- Models material memory and damage evolution
- Physically motivated by continuum damage mechanics
- Exponential approach to steady state

**When to Use:**

- Advanced material models
- Physically accurate damage evolution
- Research and detailed analysis

**Math Details:**

The transient function typically follows:

.. math::

   D(t) = D_{total} \\times \\left(1 - e^{-t/\\tau}\\right)

Where ``τ`` is a time constant related to material properties.

10. RFC_SD_TRANSIENT_23c (Recommended)
---------------------------------------

**Description:** Transient spreading from P2 to P3c (cycle closure point).

**Algorithm:**

Similar to TRANSIENT_23, but stops at the exact point P3c where the cycle closes (between P3 and P4):

.. code-block:: text

   P3c = point where min(B,C) ≥ min(A,D) condition met

   Spread damage transiently from P2 to P3c
   Stop exactly at closure, not at P4

**Visualization:**

.. code-block:: text

   Damage
     │       ╭───
     │      ╱
     │     ╱
     │    ╱
     │───┴─────────>
       P2  P3  P3c Time

   Transient response ending at cycle closure

**Why Recommended:**

1. **Physical accuracy**: Damage distribution ends exactly when cycle closes
2. **Precise timing**: Accounts for sub-sample closure position
3. **Best of both worlds**: Transient physics + exact closure timing
4. **Robust**: Works well across various load patterns

**When to Use:**

- **Default choice for most applications**
- Damage history with physical meaning desired
- Accurate timing of damage accumulation critical
- Research and production use

Implementation
==============

Enabling Damage Spreading
--------------------------

**Compile-Time:**

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_DH_SUPPORT=ON

**C Code:**

.. code-block:: c

   #include "rainflow.h"

   rfc_ctx_s ctx;
   size_t data_len = 10000;
   double data[10000];  // Your load history

   // Initialize rainflow
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Configure Wöhler curve (required for damage calculation)
   RFC_wl_init_elementary(&ctx, 1000.0, 1e7, 5.0);

   // Initialize damage history with spreading method
   RFC_dh_init(&ctx,
               RFC_SD_TRANSIENT_23c,  // Recommended method
               NULL,                  // Auto-allocate
               data_len,              // Must match data length
               false                  // Dynamic allocation
   );

   // Feed data
   RFC_feed(&ctx, data, data_len);
   RFC_finalize(&ctx, RFC_RES_REPEATED);

   // Get damage history
   const double *dh;
   size_t dh_count;
   RFC_dh_get(&ctx, &dh, &dh_count);

   // dh[i] = cumulative damage at sample i
   for (size_t i = 0; i < dh_count; i++)
   {
       printf("Sample %zu: Damage = %.6e\n", i, dh[i]);
   }

   RFC_deinit(&ctx);

**C++ Code:**

.. code-block:: cpp

   #include "rainflow.hpp"

   Rainflow::Rainflow rf;
   std::vector<double> data(10000);  // Your load history

   rf.init(100, 1.0, 0.0, 1.0);
   rf.wl_init_elementary(1000.0, 1e7, 5.0);

   // Enable damage history with spreading
   rf.dh_init(
       Rainflow::RFC_SD_TRANSIENT_23c,
       nullptr,
       data.size(),
       false
   );

   rf.feed(data.data(), data.size());
   rf.finalize(Rainflow::RFC_RES_REPEATED);

   // Get damage history
   const double* dh_data;
   size_t dh_count;
   rf.dh_get(&dh_data, &dh_count);

   std::vector<double> damage_history(dh_data, dh_data + dh_count);

**Python Example:**

.. code-block:: python

   import rfcnt
   import numpy as np

   data = np.loadtxt('measurement.txt')

   result = rfcnt.rfc(
       data,
       class_width=10.0,
       wl={'sx': 1000, 'nx': 1e7, 'k': 5},
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c  # Use spreading
   )

   damage_history = result['dh']  # Cumulative damage at each sample

Comparison of Methods
=====================

Visual Comparison
-----------------

For the same cycle B-C, different methods produce different damage distributions:

.. code-block:: text

   Signal:  A────B────C────D
                ↑    ↑
              Peak Valley

   HALF_23:           [0, D/2, D/2, 0]
   RAMP_DAMAGE_23:    [0, D/3, 2D/3, 0]
   FULL_P2:           [0, D, 0, 0]
   FULL_P3:           [0, 0, D, 0]
   TRANSIENT_23c:     [0, 0.15D, 0.60D, 0.25D]  (approx)

Damage Distribution Spread
---------------------------

.. code-block:: text

   Method              P2    P3    P4    Distribution Width
   ─────────────────────────────────────────────────────────
   HALF_23             50%   50%   0%    Narrow (2 points)
   FULL_P2             100%  0%    0%    Point (1 point)
   FULL_P3             0%    100%  0%    Point (1 point)
   RAMP_DAMAGE_23      varies      0%    Medium (P2→P3)
   RAMP_DAMAGE_24      varies varies     Wide (P2→P4)
   TRANSIENT_23c       Exponential       Medium-Wide

Computational Cost
------------------

.. code-block:: text

   Method              CPU Cost   Memory    Accuracy
   ───────────────────────────────────────────────────
   NONE                Minimal    None      N/A
   HALF_23             Low        Low       Low
   FULL_P2/P3          Low        Low       Low
   RAMP_DAMAGE_23      Low        Medium    Medium
   RAMP_AMPLITUDE_23   Medium     Medium    Medium
   RAMP_*_24           Medium     Medium    Medium
   TRANSIENT_23        High       Medium    High
   TRANSIENT_23c       High       Medium    Highest

Choosing a Spreading Method
============================

Decision Tree
-------------

.. code-block:: text

   Start
     │
     ├─ Need physically accurate damage evolution?
     │    └─ Yes → RFC_SD_TRANSIENT_23c (recommended)
     │
     ├─ Simple symmetric distribution?
     │    └─ Yes → RFC_SD_HALF_23
     │
     ├─ Damage at peak stress only?
     │    └─ Yes → RFC_SD_FULL_P2
     │
     ├─ Damage at valley/return?
     │    └─ Yes → RFC_SD_FULL_P3
     │
     ├─ Time-based linear accumulation?
     │    └─ Yes → RFC_SD_RAMP_DAMAGE_23 or _24
     │
     ├─ Amplitude-dependent accumulation?
     │    └─ Yes → RFC_SD_RAMP_AMPLITUDE_23 or _24
     │
     └─ Default → RFC_SD_TRANSIENT_23c

Recommendations by Application
-------------------------------

**General Engineering:**
   ``RFC_SD_TRANSIENT_23c`` - Best physical model

**Quick Analysis:**
   ``RFC_SD_HALF_23`` - Fast, simple

**Maximum Stress Design:**
   ``RFC_SD_FULL_P2`` - Peak-oriented

**Visualization / Plotting:**
   ``RFC_SD_TRANSIENT_23c`` or ``RFC_SD_RAMP_DAMAGE_23`` - Smooth curves

**Research / Publications:**
   ``RFC_SD_TRANSIENT_23c`` - Most defensible physically

**Performance-Critical:**
   ``RFC_SD_HALF_23`` or ``RFC_SD_FULL_P2`` - Minimal overhead

Custom Spreading via Delegates
===============================

For complete control, implement your own spreading function:

**C Code:**

.. code-block:: c

   void my_spread_damage(rfc_ctx_s *ctx,
                         rfc_value_tuple_s *from,   // P2
                         rfc_value_tuple_s *to,     // P3
                         rfc_value_tuple_s *next,   // P4
                         rfc_flags_e flags)
   {
       // Custom damage spreading logic
       double damage = /* cycle damage */;

       // Your algorithm here
       // Example: Spread based on custom material model

       // Increment damage at specific turning points
       RFC_tp_inc_damage(ctx, from->pos, damage * 0.7);
       RFC_tp_inc_damage(ctx, to->pos, damage * 0.3);
   }

   // Usage
   ctx.spread_damage_fcn = my_spread_damage;

See `delegates.rst <delegates.rst>`_ for more details on delegate functions.

Advanced Topics
===============

Damage History Resolution
--------------------------

Damage history resolution depends on input data resolution:

.. code-block:: text

   Input: 10,000 samples
   Damage History: 10,000 values (1:1 correspondence)

   Input: 1,000,000 samples
   Damage History: 1,000,000 values (8 MB memory)

For very long measurements, consider:

- Downsampling input data
- Processing in chunks
- Using compressed damage history (store only changes)

Spread Damage and Residue
--------------------------

Residue handling affects damage spreading:

.. code-block:: python

   # Residue damage is added at the end
   result = rfcnt.rfc(
       data,
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
       residue_method='repeated'
   )

   # Damage from residue appears near end of damage_history
   import matplotlib.pyplot as plt
   plt.plot(result['dh'])
   plt.xlabel('Sample')
   plt.ylabel('Cumulative Damage')
   plt.show()

Interpolation for Visualization
--------------------------------

Damage history is often piece-wise constant. Smooth for visualization:

.. code-block:: python

   from scipy.interpolate import interp1d

   # Smooth damage history
   f = interp1d(range(len(dh)), dh, kind='cubic')
   t_smooth = np.linspace(0, len(dh) - 1, len(dh) * 10)
   dh_smooth = f(t_smooth)

   plt.plot(t_smooth, dh_smooth, label='Smoothed')
   plt.plot(dh, 'o', markersize=2, label='Original')
   plt.legend()
   plt.show()

Troubleshooting
===============

Common Issues
-------------

**1. Damage history all zeros**

   .. code-block:: text

      Problem: dh array contains all zeros.

   **Causes:**

   - RFC_DH_SUPPORT not enabled at compile time
   - spread_damage method set to RFC_SD_NONE
   - Wöhler curve not initialized (no damage calculation)
   - No cycles detected

   **Solutions:**

   .. code-block:: bash

      # Ensure compiled with DH support
      cmake -S. -Bbuild -DRFC_DH_SUPPORT=ON

   .. code-block:: c

      // Use non-NONE spreading method
      RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, NULL, data_len, false);

      // Initialize Wöhler curve
      RFC_wl_init_elementary(&ctx, 1000.0, 1e7, 5.0);

**2. Damage concentrated at end**

   .. code-block:: text

      Problem: All damage appears at end of damage_history.

   **Cause:** Using RFC_SD_FULL_P3 or similar point-based method with specific data pattern.

   **Solution:** Use a spreading method like RFC_SD_TRANSIENT_23c for distributed damage.

**3. Damage history length mismatch**

   .. code-block:: text

      Error: dh_count (9999) != data_len (10000)

   **Cause:** ``data_len`` parameter in ``RFC_dh_init()`` doesn't match actual data length.

   **Solution:**

   .. code-block:: c

      size_t data_len = 10000;
      RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, NULL, data_len, false);
      RFC_feed(&ctx, data, data_len);  // Must match!

**4. Unexpected damage distribution**

   .. code-block:: text

      Problem: Damage doesn't match expectations.

   **Debugging:**

   .. code-block:: c

      // Enable debug output
      #define RFC_DEBUG_FLAGS
      ctx.debug_flags = RFC_DEBUG_DH;  // Debug damage history

      // Check each cycle's damage contribution
      // (implementation-specific)

See Also
========

- `damage_history.rst <damage_history.rst>`_ - Complete damage history documentation
- `delegates.rst <delegates.rst>`_ - Custom spreading via spread_damage_fcn delegate
- `algorithm.rst <algorithm.rst>`_ - How cycles are detected (P1, P2, P3, P4 points)
- `references.rst <references.rst>`_ - Continuum damage mechanics references

References
==========

Spread damage methods are inspired by material science and damage mechanics:

**Continuum Damage Mechanics:**

- Lemaitre, J. & Chaboche, J.L. (1990) "Mechanics of Solid Materials"
- Krajcinovic, D. (1996) "Damage Mechanics"

**Fatigue Damage Evolution:**

- Brokate, M. et al. (1996) "Rainflow counting and energy dissipation in elastoplasticity"
  Eur. J. Mech. A/Solids 15, pp. 705-737

**Transient Damage Models:**

- Hack, M. (1998) "Schaedigungsbasierte Hysteresefilter"
  Dissertation, Universität Kaiserslautern, Shaker Verlag

For implementation details:

- ``rainflow.h`` - Lines 447-460 (spread damage enum)
- ``rainflow.c`` - Function spread_damage() implementation
- ``damage_history.rst`` - Usage in damage history context
