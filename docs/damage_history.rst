================
Damage History
================

Overview
========

Damage History is a powerful feature that tracks the accumulated fatigue
damage at each point in time as the load signal is processed. This provides
insight into when and how damage accumulates during the measurement, enabling
detailed fatigue analysis and visualization.

What is Damage History?
========================

Concept
-------

When processing a time-series signal for fatigue analysis, **damage history**
records the damage increment contributed at each sample point in the input data.
Each element ``dh[i]`` holds the damage assigned to sample ``i`` by the
spreading method — most entries are zero; non-zero values appear only where
cycles are attributed.

This creates a timeline showing:

- **When** damage occurs during the measurement
- **Which parts** of the signal contribute most to total damage
- **Critical periods** where damage is concentrated

Example Visualization
---------------------

Consider a load signal with damage history::

    Load Signal:           Damage History (increments):
    ┌─────┐                |  |      |  |
    │     │                |  |      |  |
    │     │      →         |  |  |   |  |
    └─────┘                |__|__|___|__|___> Time
       Time                 (spikes at cycle-closing turning points)

The damage history shows non-zero increments at the turning points where
cycles are closed and damage is attributed. To obtain cumulative damage,
compute a prefix sum (e.g. ``np.cumsum(dh)``).

Why Use Damage History?
========================

Applications
------------

**1. Identify Critical Time Periods**

   Find when the most damaging events occur in your measurement:

   - Startup/shutdown sequences
   - Specific operational modes
   - Transient events
   - Peak load conditions

**2. Optimize Testing**

   Shorten accelerated life tests by focusing on high-damage periods:

   - Extract only damaging portions
   - Repeat critical sections
   - Skip low-damage intervals

**3. Predictive Maintenance**

   Monitor damage accumulation in real-time:

   - Estimate remaining life continuously
   - Trigger maintenance before failure
   - Optimize inspection intervals

**4. Root Cause Analysis**

   Correlate damage with operational data:

   - Temperature, speed, pressure
   - Operating modes
   - Environmental conditions

**5. Load Spectrum Editing**

   Create representative shorter test sequences:

   - Identify equivalent damage profiles
   - Compress long measurements
   - Maintain damage characteristics

Implementation
==============

Enabling Damage History
------------------------

**Compile-Time:**

The feature must be enabled during compilation:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_DH_SUPPORT=1

**C Code:**

.. code-block:: c

   #include "rainflow.h"

   rfc_ctx_s ctx;
   size_t data_len = 10000;

   // Initialize rainflow context
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Configure Wöhler curve
   RFC_wl_init_modified(&ctx, 1000.0, 1e7, 5.0, 5.0);

   // Initialize damage history
   // spread_damage: method for distributing damage to turning points
   // stream: NULL for automatic allocation
   // data_len: length of input data (must match)
   // is_static: false for dynamic allocation
   RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, NULL, data_len, false);

   // Feed data
   RFC_feed(&ctx, data, data_len);
   RFC_finalize(&ctx, RFC_RES_REPEATED);

   // Get damage history
   const double *dh;
   size_t dh_count;
   RFC_dh_get(&ctx, &dh, &dh_count);

   // Access damage increment at each time point
   // Most entries are 0; non-zero where cycles are attributed
   for (size_t i = 0; i < dh_count; i++)
   {
       if (dh[i] > 0.0)
           printf("Sample %zu: Damage increment = %.6e\n", i, dh[i]);
   }

   RFC_deinit(&ctx);

**C++ Code:**

.. code-block:: cpp

   #include "rainflow.hpp"

   Rainflow::Rainflow rf;
   size_t data_len = 10000;

   rf.init(100, 1.0, 0.0, 1.0);
   rf.wl_init_modified(1000.0, 1e7, 5.0, 5.0);

   // Enable damage history
   rf.dh_init(
       Rainflow::RFC_SD_TRANSIENT_23c,   // Spread method
       nullptr,                          // Auto-allocate
       data_len,                         // Data length
       false                             // Dynamic allocation
   );

   rf.feed(data.data(), data.size());
   rf.finalize(Rainflow::RFC_RES_REPEATED);

   // Get damage history
   const double* dh_data;
   size_t dh_count;
   rf.dh_get(&dh_data, &dh_count);

   // Copy to vector
   std::vector<double> damage_history(dh_data, dh_data + dh_count);

**Python:**

.. code-block:: python

   import rfcnt
   import numpy as np

   data = np.loadtxt('measurement.txt')

   result = rfcnt.rfc(
       data,
       class_width=10.0,
       wl={'sx': 1000, 'nx': 1e7, 'k': 5, 'k2': 5},
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c  # Enable DH
   )

   # Access damage history
   damage_history = result['dh']

   # Plot
   import matplotlib.pyplot as plt
   plt.plot(damage_history)
   plt.xlabel('Sample Index')
   plt.ylabel('Accumulated Damage')
   plt.title('Damage History')
   plt.show()

Damage Spreading Methods
=========================

The ``spread_damage`` parameter controls how damage from closed cycles is
distributed to the turning points. See `spread_damage.rst <spread_damage.rst>`_
for detailed explanation of each method.

Available Methods
-----------------

.. code-block:: c

   enum rfc_sd_method_e
   {
       RFC_SD_NONE             = -1, // No damage history
       RFC_SD_HALF_23          =  0, // Equally split between P2 and P3
       RFC_SD_RAMP_AMPLITUDE_23,     // Ramp by amplitude over P2 to P3
       RFC_SD_RAMP_DAMAGE_23,        // Ramp evenly over P2 to P3
       RFC_SD_RAMP_AMPLITUDE_24,     // Ramp by amplitude over P2 to P4
       RFC_SD_RAMP_DAMAGE_24,        // Ramp evenly over P2 to P4
       RFC_SD_FULL_P2,               // Full damage at P2
       RFC_SD_FULL_P3,               // Full damage at P3
       RFC_SD_TRANSIENT_23,          // Transient over P2 to P3
       RFC_SD_TRANSIENT_23c,         // Transient over P2 to P4, until cycle closes (recommended)
       RFC_SD_COUNT                  // Number of options
   };

**Recommended:** ``RFC_SD_TRANSIENT_23c`` provides the most physically
realistic damage distribution.

Storage Modes
=============

Compressed vs. Uncompressed
----------------------------

**Uncompressed Mode (Default)**:

- Stores damage for every input sample
- Array length = input data length
- Direct 1:1 correspondence with input
- Higher memory usage

**Compressed Mode**:

- Stores damage only at turning points
- Smaller memory footprint
- Must correlate with turning point positions
- More complex to visualize

Uncompressed mode is recommended for most applications as it simplifies
visualization and analysis.

Memory Requirements
-------------------

For a measurement with ``N`` samples:

.. code-block:: text

   Uncompressed: N × sizeof(double) = N × 8 bytes

   Example: 1,000,000 samples → 8 MB
            10,000,000 samples → 80 MB

For embedded systems or very long measurements, consider:

- Processing in chunks
- Using compressed mode
- Storing only critical portions

Accessing Damage History
=========================

Getting the Data
----------------

The damage history array contains per-sample damage increments:

.. code-block:: c

   const double *dh;
   size_t dh_count;

   RFC_dh_get(&ctx, &dh, &dh_count);

   // dh[i] = damage increment attributed to sample i
   // Most entries are 0; non-zero at cycle-closing turning points
   // sum(dh[0..dh_count-1]) == total damage

**Properties:**

- **Non-negative**: ``dh[i] >= 0``
- **Sparse**: most entries are zero
- **Sums to total damage**: ``sum(dh) == total_damage``
- To obtain cumulative damage over time, apply a prefix sum

Cumulative Damage and Damage Rate
----------------------------------

``dh`` contains damage increments. To get cumulative damage over time:

.. code-block:: python

   import numpy as np

   # dh contains increments — cumulate to get damage over time
   cumulative_damage = np.cumsum(dh)

   # dh itself is already the damage rate per sample
   # Find samples with highest damage contribution
   critical_indices = np.argsort(dh)[-100:]

   print(f"Most damaging samples: {critical_indices}")

Visualization
=============

Basic Plot
----------

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np

   # Assume damage_history from rfcnt.rfc()
   time = np.arange(len(damage_history))

   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

   # Plot 1: Damage increments
   ax1.plot(time, damage_history, 'b-', linewidth=1)
   ax1.set_xlabel('Sample')
   ax1.set_ylabel('Damage Increment')
   ax1.set_title('Damage History (increments per sample)')
   ax1.grid(True)

   # Plot 2: Cumulative damage over time
   cumulative = np.cumsum(damage_history)
   ax2.plot(time, cumulative, 'r-', linewidth=1)
   ax2.set_xlabel('Sample')
   ax2.set_ylabel('Cumulative Damage')
   ax2.set_title('Cumulative Damage')
   ax2.grid(True)

   plt.tight_layout()
   plt.show()

Correlation with Load Signal
-----------------------------

.. code-block:: python

   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

   # Original signal
   ax1.plot(data, 'k-', linewidth=0.5, label='Load Signal')
   ax1.set_ylabel('Load')
   ax1.legend()
   ax1.grid(True)

   # Damage history
   ax2.plot(damage_history, 'r-', linewidth=1, label='Damage History')
   ax2.set_xlabel('Sample')
   ax2.set_ylabel('Cumulative Damage')
   ax2.legend()
   ax2.grid(True)

   plt.tight_layout()
   plt.show()

Highlighting Critical Periods
------------------------------

.. code-block:: python

   # Find regions with high damage (dh already contains increments)
   threshold = np.percentile(dh[dh > 0], 75)  # Top 25% of non-zero increments

   plt.figure(figsize=(14, 6))
   plt.plot(data, 'k-', linewidth=0.5, alpha=0.5, label='Load Signal')

   # Highlight critical regions
   critical = dh > threshold
   plt.fill_between(
       range(len(data)),
       data.min(), data.max(),
       where=critical,
       color='red',
       alpha=0.3,
       label='Critical Damage Periods'
   )

   plt.xlabel('Sample')
   plt.ylabel('Load')
   plt.title('Load Signal with Critical Damage Periods Highlighted')
   plt.legend()
   plt.grid(True)
   plt.show()

Practical Applications
======================

Example 1: Test Optimization
-----------------------------

Identify the most damaging 10% of a measurement:

.. code-block:: python

   import rfcnt
   import numpy as np

   # Original measurement (very long)
   data = np.loadtxt('long_measurement.txt')  # 1,000,000 samples

   # Compute damage history
   result = rfcnt.rfc(
       data,
       class_width=10.0,
       wl={'sx': 1000, 'nx': 1e7, 'k': 5},
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c
   )

   dh = result['dh']
   total_damage = dh.sum()
   cumulative = np.cumsum(dh)

   # Find samples contributing to the middle 90% of cumulative damage
   critical_samples = np.where(
       (cumulative >= 0.05 * total_damage) &  # After 5% accumulated
       (cumulative <= 0.95 * total_damage)     # Before 95% accumulated
   )[0]

   # Extract critical portion
   optimized_test = data[critical_samples]

   print(f"Original length: {len(data)} samples")
   print(f"Optimized length: {len(optimized_test)} samples")
   print(f"Reduction: {100*(1-len(optimized_test)/len(data)):.1f}%")
   print(f"Damage retained: 90%")

Example 2: Remaining Life Estimation
-------------------------------------

.. code-block:: python

   def estimate_remaining_life(damage_history, total_cycles=1):
       """
       Estimate remaining life assuming Miner's rule.

       Parameters:
       -----------
       damage_history : array
           Cumulative damage over time
       total_cycles : int
           Number of times this load sequence will repeat

       Returns:
       --------
       remaining_cycles : float
           Estimated remaining repetitions before failure
       """
       total_damage = damage_history.sum()
       damage_per_cycle = total_damage / total_cycles

       # Failure at D = 1.0 (Miner's rule)
       remaining_damage = 1.0 - total_damage
       remaining_cycles = remaining_damage / damage_per_cycle

       return remaining_cycles

   # Usage
   remaining = estimate_remaining_life(result['dh'])
   print(f"Estimated remaining life: {remaining:.0f} cycles")

Example 3: Real-Time Monitoring (C++)
--------------------------------------

Stream data in chunks using ``feed()``, reading the scalar damage accumulator
after each chunk. No damage history array is allocated — O(1) memory overhead.
The ``ctx_get().damage`` field is updated whenever a cycle closes, so it
reflects the damage from all turning points processed so far.

.. code-block:: cpp

   #include "rainflow.hpp"
   #include <iostream>
   #include <vector>

   // User-supplied: fill buf with next chunk, return actual count (0 = done)
   size_t get_next_chunk( double* buf, size_t capacity );

   int main()
   {
       const double DAMAGE_LIMIT     = 1.0;
       const double WARN_THRESHOLD   = 0.8;
       const size_t CHUNK_SIZE       = 256;

       Rainflow rf;

       if( !rf.init( 100, 1.0, 0.0, 1.0 ) )
           return 1;

       // Wöhler curve required for damage calculation
       rf.wl_init_modified( 1000.0, 1e7, 5.0, 5.0 );

       // No dh_init() — damage is tracked as a scalar, not per sample

       std::vector<double> chunk( CHUNK_SIZE );
       size_t total_samples = 0;
       bool alarm = false;

       while( !alarm )
       {
           size_t n = get_next_chunk( chunk.data(), CHUNK_SIZE );
           if( n == 0 ) break;

           if( !rf.feed( chunk.data(), n ) )
           {
               std::cerr << "feed() error: " << rf.error_get() << std::endl;
               return 1;
           }

           total_samples += n;

           // damage accumulates as cycles close at turning points
           double current_damage = rf.ctx_get().damage;
           double remaining_pct  = 100.0 * ( 1.0 - current_damage / DAMAGE_LIMIT );

           if( current_damage >= DAMAGE_LIMIT )
           {
               std::cout << "Sample " << total_samples
                         << ": FAILURE PREDICTED (D=" << current_damage << ")\n";
               alarm = true;
           }
           else if( current_damage >= WARN_THRESHOLD * DAMAGE_LIMIT )
           {
               std::cout << "Sample " << total_samples
                         << ": WARNING — " << remaining_pct << "% life remaining\n";
           }
       }

       // Finalize: process residue (unclosed cycles at end of stream)
       rf.finalize( Rainflow::RFC_RES_REPEATED );

       double total_damage;
       rf.damage( &total_damage );
       std::cout << "Final damage (incl. residue): " << total_damage << "\n";

       return 0;
   }

Performance Considerations
==========================

Computational Cost
------------------

Damage history adds minimal computational overhead:

- **Memory**: O(N) additional storage
- **Time**: O(1) per cycle (just accumulation)
- **Impact**: < 5% performance degradation

Memory Management
-----------------

For large datasets, consider:

**1. Chunked Processing:**

   .. warning::

      Each ``rfcnt.rfc()`` call is independent: cycles that span a chunk
      boundary are never closed, so their damage is lost. This underestimates
      total damage. For accurate results use the C++ ``feed()`` streaming API
      (see Example 3 above) or process the full signal in one call.

   .. code-block:: python

      chunk_size = 100000
      total_damage = 0.0

      for chunk_start in range(0, len(data), chunk_size):
          chunk = data[chunk_start:chunk_start + chunk_size]
          result = rfcnt.rfc(chunk, class_width=10.0, ...)
          total_damage += result['damage']

**2. Downsampling DH:**

   .. code-block:: python

      # Sum increments within each block of 100 samples — no damage is lost
      dh_blocked = np.add.reduceat(damage_history, np.arange(0, len(damage_history), 100))

**3. Selective Storage:**

   .. code-block:: python

      # Store only samples with non-negligible damage increments
      threshold = np.max(damage_history) * 1e-6  # relative to peak increment
      keep = damage_history > threshold
      sparse_indices = np.where(keep)[0]
      sparse_dh = damage_history[keep]  # (index, increment) pairs

Limitations
===========

Current Limitations
-------------------

1. **Sequential Processing Only**

   Damage history must be computed in order - cannot be parallelized across
   the time axis.

2. **Memory Scaling**

   Memory usage scales linearly with input length in uncompressed mode.

3. **Post-Processing Only**

   The full damage history is only available after all data is processed
   (use streaming for real-time monitoring).

4. **Residue Handling**

   Damage from residue (if any) is added at the end, not distributed
   throughout the signal.

Troubleshooting
===============

Common Issues
-------------

**1. DH array length doesn't match input**

   .. code-block:: text

      Error: dh_count (9999) != data_len (10000)

   **Solution**: Ensure ``data_len`` parameter in ``RFC_dh_init()`` exactly
   matches the actual data length you'll feed.

**2. All damage at the end**

   **Cause**: Wrong spread method or residue-only damage.

   **Solution**: Use ``RFC_SD_TRANSIENT_23c`` and ensure cycles are being
   closed during processing, not just in residue.

**3. Out of memory**

   **Cause**: Very long measurement.

   **Solution**: Process in chunks or use compressed mode.

See Also
========

- `spread_damage.rst <spread_damage.rst>`_ - Damage spreading methods explained
- `examples.rst <examples.rst>`_ - More code examples
- `features.rst <features.rst>`_ - Feature compilation options
- `residue_methods.rst <residue_methods.rst>`_ - Residue handling