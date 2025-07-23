===============
Turning Points
===============

Overview
========

Turning points are the local maxima and minima in a load signal - the peaks and valleys that define fatigue cycles. This document explains how the rainflow library handles turning point storage, including internal vs external storage options, and how to work with turning point data.

What are Turning Points?
=========================

Definition
----------

A **turning point** is a point in the load-time history where the signal changes direction:

- **Peak (local maximum)**: Load increases then decreases
- **Valley (local minimum)**: Load decreases then increases

Example Signal::

    Load
     │
     │       *B      *D
     │      / \    / \
     │     /   \  /   \       B, D = Peaks (turning points)
     │    /     \/     \
     │   /      *C      \ /   A, C, E = Valleys (turning points)
     │\ /                *E
     │ *A
     └─────────────────────────────> Time
                           (E is above A level)

The signal has 5 turning points: A (valley), B (peak), C (valley), D (peak), E (valley).

Importance in Rainflow Counting
--------------------------------

Rainflow counting operates exclusively on turning points, not the raw signal:

1. **Dimensionality reduction**: 10,000 data points might yield only 100 turning points
2. **Cycle detection**: The 4-point algorithm examines sequences of 4 turning points
3. **Storage efficiency**: Only turning points need to be stored, not all raw data

**Example:**

.. code-block:: text

   Raw data:     [0, 5, 10, 8, 10, 9, 10, 5, 0]  (9 points)
   Turning pts:  [0, 10, 5, 0]                   (4 points)

The rainflow algorithm sees only the 4 turning points.

RFC_TP_SUPPORT Feature
=======================

Enabling Turning Point Storage
-------------------------------

Compile with:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_TP_SUPPORT=ON

Or:

.. code-block:: c

   #define RFC_TP_SUPPORT 1

**Effect:** When enabled, the library stores all turning points internally and provides API functions to access them.

When Disabled (Default in RFC_MINIMAL)
---------------------------------------

Without ``RFC_TP_SUPPORT``:

- Turning points are detected but immediately processed
- Only the current working set (residue) is kept in memory
- Cannot retrieve full turning point history
- Lower memory usage

Storage Modes
=============

The library supports two storage modes for turning points:

1. Internal Storage
-------------------

**Description:** Library allocates and manages a buffer for turning points internally.

**How it Works:**

- During ``RFC_init()``, a buffer is allocated (``ctx->tp``)
- As turning points are detected, they're stored in this buffer
- Buffer automatically resizes if needed (up to capacity limit)

**Advantages:**

- Simple to use - no external management needed
- Direct access to turning point array
- Fast random access

**Disadvantages:**

- Memory usage scales with number of turning points
- Limited by available RAM (problematic for embedded systems)
- All turning points must fit in memory simultaneously

**Example:**

.. code-block:: c

   #include "rainflow.h"

   rfc_ctx_s ctx;
   double data[] = {0, 10, 0, 20, 0, 30, 0};

   // Initialize with internal TP storage
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Feed data - turning points stored internally
   RFC_feed(&ctx, data, 7);

   // Access turning points
   printf("Stored %zu turning points\n", ctx.tp_cnt);

   for (size_t i = 0; i < ctx.tp_cnt; i++)
   {
       printf("TP %zu: value=%.1f, class=%u, pos=%zu\n",
              i,
              ctx.tp[i].value,
              ctx.tp[i].cls,
              ctx.tp[i].pos);
   }

   RFC_deinit(&ctx);

**Output:**

.. code-block:: text

   Stored 7 turning points
   TP 0: value=0.0, class=50, pos=1
   TP 1: value=10.0, class=60, pos=2
   TP 2: value=0.0, class=50, pos=3
   TP 3: value=20.0, class=70, pos=4
   TP 4: value=0.0, class=50, pos=5
   TP 5: value=30.0, class=80, pos=6
   TP 6: value=0.0, class=50, pos=7

2. External Storage
-------------------

**Description:** You provide custom storage and access functions via delegates.

**How it Works:**

- Set ``tp_set_fcn`` and ``tp_get_fcn`` delegates
- Library calls your functions to read/write turning points
- You control where and how data is stored

**Advantages:**

- Unlimited storage capacity (database, file, external RAM)
- Flexible: can use circular buffers, compression, etc.
- Suitable for embedded systems with limited RAM

**Disadvantages:**

- More complex to implement
- Potentially slower (I/O overhead)
- Requires careful management of storage lifetime

**Example:**

.. code-block:: c

   #include "rainflow.h"
   #include <stdio.h>

   // External storage in a file
   FILE *tp_file = NULL;
   rfc_value_tuple_s tp_cache[10];  // Small cache

   bool tp_set_external(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
   {
       // Write to file
       fseek(tp_file, (tp_pos - 1) * sizeof(*tp), SEEK_SET);
       fwrite(tp, sizeof(*tp), 1, tp_file);

       // Also cache
       if (tp_pos <= 10)
       {
           tp_cache[tp_pos - 1] = *tp;
       }

       return true;
   }

   bool tp_get_external(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s **tp)
   {
       // Read from cache or file
       if (tp_pos <= 10)
       {
           *tp = &tp_cache[tp_pos - 1];
       }
       else
       {
           static rfc_value_tuple_s temp;
           fseek(tp_file, (tp_pos - 1) * sizeof(temp), SEEK_SET);
           fread(&temp, sizeof(temp), 1, tp_file);
           *tp = &temp;
       }

       return true;
   }

   int main(void)
   {
       rfc_ctx_s ctx;
       double data[] = {0, 10, 0, 20, 0, 30, 0};

       // Open external storage
       tp_file = fopen("turning_points.bin", "wb+");

       // Initialize rainflow
       RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

       // Set delegates for external storage
       ctx.tp_set_fcn = tp_set_external;
       ctx.tp_get_fcn = tp_get_external;

       // Set tp to NULL to indicate external storage
       ctx.tp = NULL;

       // Feed data - turning points stored externally
       RFC_feed(&ctx, data, 7);

       printf("Stored %zu turning points externally\n", ctx.tp_cnt);

       // Cleanup
       RFC_deinit(&ctx);
       fclose(tp_file);

       return 0;
   }

See `delegates.rst <delegates.rst>`_ for more details on delegate functions.

Turning Point Data Structure
=============================

Each turning point is represented by ``rfc_value_tuple_s``:

.. code-block:: c

   struct rfc_value_tuple
   {
       rfc_value_t  value;    // Load value (double)
       unsigned     cls;      // Discretized class number (base 0)
       size_t       pos;      // Position in input stream (1-based)

   #if RFC_TP_SUPPORT
       size_t       adj_pos;  // Position of adjacent paired turning point (1-based)
       size_t       tp_pos;   // Position in tp storage (1-based); 0 in tp array itself
       rfc_value_t  avrg;     // Average value of the two paired turning points
   #if RFC_DH_SUPPORT
       double       damage;   // Damage accumulated to this turning point
   #endif
   #endif
   };

**Fields:**

- **value**: The actual load value at the turning point
- **cls**: The discretized class (bin) number based on ``class_width`` (0-based)
- **pos**: Position in the original input data stream (1-based indexing)
- **adj_pos**: (if ``RFC_TP_SUPPORT``) Position of the adjacent paired turning point
- **tp_pos**: (if ``RFC_TP_SUPPORT``) Position of this point in the tp storage; always 0 inside the ``ctx.tp`` array itself
- **avrg**: (if ``RFC_TP_SUPPORT``) Average value of the two paired turning points
- **damage**: (if ``RFC_TP_SUPPORT`` && ``RFC_DH_SUPPORT``) Damage accumulated at this point

**Example:**

.. code-block:: c

   rfc_value_tuple_s tp = ctx.tp[3];

   printf("Value: %.2f\n", tp.value);      // 20.00
   printf("Class: %u\n", tp.cls);          // 70
   printf("Position: %zu\n", tp.pos);      // 4
   printf("Damage: %.3e\n", tp.damage);    // 0.000123

API Functions
=============

When ``RFC_TP_SUPPORT`` is enabled, these functions are available:

1. RFC_tp_init()
----------------

Initialize turning point storage with a pre-allocated buffer.

**Signature:**

.. code-block:: c

   bool RFC_tp_init(
       rfc_ctx_s         *rfc_ctx,
       rfc_value_tuple_s *tp,        // Buffer (or NULL for dynamic allocation)
       size_t             tp_cap,    // Initial capacity
       bool               is_static  // true if buffer is statically allocated
   );

**Example:**

.. code-block:: c

   // Static buffer
   static rfc_value_tuple_s tp_buf[1000];
   RFC_tp_init(&ctx, tp_buf, 1000, true);

2. RFC_tp_prune()
-----------------

Remove turning points from storage to free memory.

See `tp_prune.rst <tp_prune.rst>`_ for detailed explanation of this complex function.

3. RFC_tp_refeed()
------------------

Re-process stored turning points through rainflow counting, optionally with
changed hysteresis or class parameters.

**Signature:**

.. code-block:: c

   bool RFC_tp_refeed(
       rfc_ctx_s               *rfc_ctx,
       rfc_value_t              new_hysteresis,   // New hysteresis (0 to keep current)
       const rfc_class_param_s *new_class_param   // New class params (NULL to keep current)
   );

**Purpose:** Re-run rainflow counting on already-stored turning points without
feeding raw data again.

**Use Cases:**

- Reprocessing with different hysteresis or class parameters
- Testing different residue methods via ``RFC_finalize()``

**Example:**

.. code-block:: c

   // Initial processing
   RFC_feed(&ctx, data, data_len);
   RFC_finalize(&ctx, RFC_RES_IGNORE);

   // Reprocess with the same parameters but a different residue method
   RFC_tp_refeed(&ctx, 0.0, NULL);
   RFC_finalize(&ctx, RFC_RES_REPEATED);

4. RFC_tp_clear()
-----------------

Clear the turning point buffer (resets ``ctx.tp_cnt`` to 0; does not free memory).

**Signature:**

.. code-block:: c

   bool RFC_tp_clear( rfc_ctx_s *rfc_ctx );

Accessing Individual Turning Points
-------------------------------------

There is no public ``RFC_tp_get()`` or ``RFC_tp_set()`` function. Access turning
points directly or via the delegate fields:

**Internal storage** — direct array access (0-based index):

.. code-block:: c

   for (size_t i = 0; i < ctx.tp_cnt; i++)
   {
       printf("TP %zu: value=%.1f, class=%u, pos=%zu\n",
              i, ctx.tp[i].value, ctx.tp[i].cls, ctx.tp[i].pos);
   }

**External storage** — call the delegates (1-based position):

.. code-block:: c

   rfc_value_tuple_s *tp;
   if (ctx.tp_get_fcn && ctx.tp_get_fcn(&ctx, i + 1, &tp))
   {
       printf("TP %zu: value=%.1f\n", i + 1, tp->value);
   }

To write a turning point, use ``ctx.tp_set_fcn`` or direct array write.
To increment damage at a turning point, use ``ctx.tp_inc_damage_fcn``.

Accessing Turning Point Data
=============================

Direct Access (Internal Storage)
---------------------------------

When using internal storage, you can directly access the array:

.. code-block:: c

   rfc_ctx_s ctx;
   // ... initialize and feed data ...

   printf("Total turning points: %zu\n", ctx.tp_cnt);

   // Direct array access
   for (size_t i = 0; i < ctx.tp_cnt; i++)
   {
       rfc_value_tuple_s *tp = &ctx.tp[i];
       printf("TP %zu: value=%.2f\n", i, tp->value);
   }

**Caution:** Modifying ``ctx.tp`` directly can corrupt internal state. Use API functions when possible.

Iterator-Style Access
---------------------

For portability (works with both internal and external storage), use ``RFC_tp_get()``:

.. code-block:: c

   for (size_t i = 1; i <= ctx.tp_cnt; i++)  // Note: 1-based
   {
       rfc_value_tuple_s *tp;
       if (RFC_tp_get(&ctx, i, &tp))
       {
           printf("TP %zu: value=%.2f\n", i, tp->value);
       }
   }

Exporting Turning Points
-------------------------

Export to CSV:

.. code-block:: c

   void export_turning_points_csv(rfc_ctx_s *ctx, const char *filename)
   {
       FILE *fp = fopen(filename, "w");
       if (!fp) return;

       fprintf(fp, "Index,Value,Class,Position,Damage\n");

       for (size_t i = 1; i <= ctx->tp_cnt; i++)
       {
           rfc_value_tuple_s *tp;
           if (RFC_tp_get(ctx, i, &tp))
           {
               fprintf(fp, "%zu,%.6f,%u,%zu,%.6e\n",
                      i, tp->value, tp->cls, tp->pos,
                      #if RFC_TP_SUPPORT && RFC_DH_SUPPORT
                      tp->damage
                      #else
                      0.0
                      #endif
               );
           }
       }

       fclose(fp);
   }

Practical Applications
======================

1. Signal Compression
---------------------

Turning points provide significant data compression:

**Example:**

.. code-block:: c

   size_t raw_data_count = 100000;
   size_t tp_count = ctx.tp_cnt;  // Might be 500

   double compression_ratio = (double)raw_data_count / tp_count;
   printf("Compression: %.1fx (100k points → %zu TPs)\n",
          compression_ratio, tp_count);

**Output:**

.. code-block:: text

   Compression: 200.0x (100k points → 500 TPs)

2. Cycle Visualization
----------------------

Plot cycles from turning points:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np

   # Load turning points from CSV export
   tps = np.loadtxt('turning_points.csv', delimiter=',', skiprows=1)

   plt.figure(figsize=(12, 6))
   plt.plot(tps[:, 3], tps[:, 1], 'o-', markersize=4)
   plt.xlabel('Position in Input Stream')
   plt.ylabel('Load Value')
   plt.title('Turning Points from Rainflow Analysis')
   plt.grid(True)
   plt.show()

3. Pre-Filtering
----------------

Apply custom filters before rainflow counting:

.. code-block:: c

   void custom_tp_filter(rfc_ctx_s *ctx)
   {
       // Example: Remove small turning points
       double threshold = 5.0;

       size_t write_idx = 0;
       for (size_t read_idx = 0; read_idx < ctx->tp_cnt; read_idx++)
       {
           rfc_value_tuple_s *tp = &ctx->tp[read_idx];

           // Keep only significant turning points
           if (read_idx == 0 ||
               read_idx == ctx->tp_cnt - 1 ||
               fabs(tp->value - ctx->tp[read_idx - 1].value) > threshold)
           {
               if (write_idx != read_idx)
               {
                   ctx->tp[write_idx] = *tp;
               }
               write_idx++;
           }
       }

       ctx->tp_cnt = write_idx;
       printf("Filtered: %zu → %zu turning points\n",
              ctx->tp_cnt, write_idx);
   }

4. Incremental Processing
--------------------------

Process turning points in chunks:

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Process data in chunks
   for (int chunk = 0; chunk < NUM_CHUNKS; chunk++)
   {
       double *chunk_data = load_chunk(chunk);
       size_t chunk_len = get_chunk_len(chunk);

       size_t tp_before = ctx.tp_cnt;
       RFC_feed(&ctx, chunk_data, chunk_len);
       size_t tp_after = ctx.tp_cnt;

       printf("Chunk %d: added %zu turning points\n",
              chunk, tp_after - tp_before);

       free(chunk_data);
   }

   RFC_finalize(&ctx, RFC_RES_REPEATED);

Memory Considerations
=====================

Memory Usage
------------

For internal storage:

.. code-block:: text

   Memory = tp_cap × sizeof(rfc_value_tuple_s)

   sizeof depends on compile flags (64-bit platform):
     Without RFC_TP_SUPPORT: value(8) + cls(4) + pad(4) + pos(8) ≈ 24 bytes
     With RFC_TP_SUPPORT:    + adj_pos(8) + tp_pos(8) + avrg(8) ≈ 48 bytes
     + RFC_DH_SUPPORT:       + damage(8)                         ≈ 56 bytes

   Example (RFC_TP_SUPPORT + RFC_DH_SUPPORT):
     10,000 TPs  ≈ 560 KB
     100,000 TPs ≈ 5.6 MB

Capacity Management
-------------------

The library automatically resizes the TP buffer, but you can set initial capacity:

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Pre-allocate for expected number of TPs
   ctx.tp_cap = 10000;  // Expect ~10k turning points
   ctx.tp = malloc(ctx.tp_cap * sizeof(rfc_value_tuple_s));

Reducing Memory Usage
---------------------

**1. Use RFC_tp_prune()**

   Periodically remove processed turning points:

   .. code-block:: c

      if (ctx.tp_cnt > 5000)
      {
          RFC_tp_prune(&ctx, 0, RFC_FLAGS_DEFAULT);  // Keep only residue
      }

**2. External Storage**

   Stream turning points to disk/database instead of RAM.

**3. Disable TP Storage**

   If you don't need turning point history:

   .. code-block:: bash

      cmake -S. -Bbuild -DRFC_TP_SUPPORT=OFF

   Memory savings: ~30 bytes per turning point.

Performance Considerations
==========================

Access Patterns
---------------

**Internal Storage:**

- **Random access**: O(1) - direct array indexing
- **Sequential access**: O(n) - cache-friendly
- **Best for**: Fast processing, abundant RAM

**External Storage:**

- **Random access**: O(1) to O(log n) - depends on backing store
- **Sequential access**: O(n) - may have I/O overhead
- **Best for**: Limited RAM, very large datasets

Optimization Tips
-----------------

**1. Pre-allocate Capacity**

   Avoid repeated reallocations:

   .. code-block:: c

      // Bad: Default small capacity, many reallocations
      RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

      // Good: Pre-allocate expected capacity
      RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
      ctx.tp_cap = expected_tp_count;
      ctx.tp = realloc(ctx.tp, ctx.tp_cap * sizeof(*ctx.tp));

**2. Batch External Writes**

   When using external storage, buffer writes:

   .. code-block:: c

      #define BUFFER_SIZE 100
      rfc_value_tuple_s buffer[BUFFER_SIZE];
      size_t buffer_count = 0;

      bool tp_set_buffered(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
      {
          buffer[buffer_count++] = *tp;

          if (buffer_count >= BUFFER_SIZE)
          {
              // Flush to external storage
              write_to_storage(buffer, buffer_count);
              buffer_count = 0;
          }

          return true;
      }

**3. Compress Turning Points**

   For very large datasets, consider compression:

   - Store deltas instead of absolute values
   - Quantize values more aggressively
   - Use run-length encoding for repeated classes

C++ Integration
===============

Using STL Containers
--------------------

Store turning points in C++ containers:

.. code-block:: cpp

   #include "rainflow.hpp"
   #include <vector>

   // External storage using std::vector
   std::vector<rfc_value_tuple_s> tp_storage;

   bool tp_set_vector(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
   {
       if (tp_pos > tp_storage.size())
       {
           tp_storage.resize(tp_pos);
       }
       tp_storage[tp_pos - 1] = *tp;
       return true;
   }

   bool tp_get_vector(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s **tp)
   {
       if (tp_pos == 0 || tp_pos > tp_storage.size())
           return false;

       *tp = &tp_storage[tp_pos - 1];
       return true;
   }

   // Usage
   Rainflow rf;
   rf.init(100, 1.0);

   rf.ctx_get().tp_set_fcn = tp_set_vector;
   rf.ctx_get().tp_get_fcn = tp_get_vector;
   rf.ctx_get().tp = nullptr;  // Signal external storage

RAII Wrapper
------------

.. code-block:: cpp

   class TurningPointManager
   {
   public:
       TurningPointManager(size_t capacity)
       {
           storage_.reserve(capacity);
       }

       bool set(size_t pos, const rfc_value_tuple_s *tp)
       {
           if (pos == 0) return false;
           if (pos > storage_.size())
               storage_.resize(pos);
           storage_[pos - 1] = *tp;
           return true;
       }

       bool get(size_t pos, rfc_value_tuple_s **tp)
       {
           if (pos == 0 || pos > storage_.size())
               return false;
           *tp = &storage_[pos - 1];
           return true;
       }

       const std::vector<rfc_value_tuple_s>& data() const
       {
           return storage_;
       }

   private:
       std::vector<rfc_value_tuple_s> storage_;
   };

Troubleshooting
===============

Common Issues
-------------

**1. tp_cnt is 0 after RFC_feed()**

   .. code-block:: text

      Problem: No turning points detected.

   **Solutions:**

   - Check that RFC_TP_SUPPORT is enabled
   - Verify input data has variations (not all constant values)
   - Check hysteresis settings (if enabled)

**2. Segmentation fault accessing ctx.tp**

   .. code-block:: text

      Problem: Crash when reading turning point array.

   **Solutions:**

   - Ensure RFC_TP_SUPPORT was enabled at compile time
   - Check that ``ctx.tp`` is not NULL
   - Use ``RFC_tp_get()`` instead of direct array access
   - Verify index is within bounds (0 to ctx.tp_cnt - 1)

**3. Out of memory**

   .. code-block:: text

      Problem: malloc() fails for large TP count.

   **Solutions:**

   - Use external storage with delegates
   - Enable RFC_tp_prune() periodic cleanup
   - Increase available heap size
   - Process data in smaller chunks

**4. Wrong tp_pos in delegates**

   .. code-block:: text

      Problem: tp_pos seems incorrect in set/get functions.

   **Remember:** tp_pos is **1-based**, not 0-based!

   .. code-block:: c

      // Correct
      array[tp_pos - 1] = *tp;  // Subtract 1 for array index

      // Wrong
      array[tp_pos] = *tp;  // Off-by-one error

**5. Turning points missing**

   .. code-block:: text

      Problem: Expected more turning points.

   **Causes:**

   - Hysteresis filtering is active (RFC_HYSTERESIS_FILTER)
   - Class width too coarse - adjacent points quantize to same class
   - Peak-valley filtering removed insignificant points

   **Solutions:**

   - Reduce ``hysteresis`` parameter
   - Reduce ``class_width`` for finer resolution
   - Check ``RFC_FLAGS_*`` settings

See Also
========

- `delegates.rst <delegates.rst>`_ - Delegate functions for external storage
- `tp_prune.rst <tp_prune.rst>`_ - Memory management with RFC_tp_prune()
- `damage_history.rst <damage_history.rst>`_ - Damage tracking at turning points
- `algorithm.rst <algorithm.rst>`_ - How turning points are used in counting
- `cpp_wrapper.rst <cpp_wrapper.rst>`_ - C++ integration examples

References
==========

Turning point detection is fundamental to rainflow counting:

- **ASTM E 1049** - Section 6.3: "Identification of Turning Points"
- **DIN 45667** - Peak-valley extraction procedures
- **Downing & Socie (1982)** - "Simple rainflow counting algorithms"

For implementation details, see:

- ``rainflow.h`` - Lines 591-595 (TP delegate typedefs)
- ``rainflow.c`` - Functions RFC_tp_get(), RFC_tp_set(), RFC_tp_refeed()
- ``rfc_test.c`` - Turning point testing examples
