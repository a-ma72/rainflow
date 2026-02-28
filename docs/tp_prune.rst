=========
TP Prune
=========

Overview
========

``RFC_tp_prune()`` is a memory management function that removes turning points from storage to prevent memory exhaustion during long measurements. This is essential for real-time monitoring systems, embedded applications, and processing very large datasets where storing all turning points would exceed available memory.

What is TP Pruning?
===================

The Problem
-----------

When processing a continuous data stream with turning point storage enabled (``RFC_TP_SUPPORT``), the turning point buffer grows indefinitely:

.. code-block:: text

   Time:      t=0 ────────────> t=∞
   Data:      Continuous stream (millions of samples)
   TPs:       [TP1, TP2, TP3, ..., TP_1000000, ...]
   Memory:    Growing... eventually exhausts RAM!

For a measurement with 1 million turning points:

.. code-block:: text

   Memory = 1,000,000 × sizeof(rfc_value_tuple_s)
          = 1,000,000 × 28 bytes
          = 28 MB

This is unacceptable for embedded systems or real-time monitoring.

The Solution: Pruning
---------------------

**Pruning** removes "old" turning points that are no longer needed for rainflow counting, keeping only:

1. **Residue** (unclosed cycles) - always needed
2. **Recent history** (optional) - for context or reprocessing

.. code-block:: text

   Before Prune:  [TP1, TP2, ..., TP1000, residue]
                   ↓ Prune to 100 points
   After Prune:   [TP900, ..., TP1000, residue]

   Memory freed: 900 × 28 bytes = 25.2 KB

Function Signature
==================

.. code-block:: c

   bool RFC_tp_prune(
       void *ctx,           // Rainflow context
       size_t limit,        // Target number of TPs after pruning
       rfc_flags_e flags    // Pruning behavior flags
   );

**Parameters:**

- ``ctx``: Pointer to ``rfc_ctx_s`` structure
- ``limit``: Desired number of turning points to keep
- ``flags``: Control flags for pruning behavior

**Returns:** ``true`` on success, ``false`` on error

**Requires:** ``RFC_TP_SUPPORT`` enabled at compile time

Pruning Modes
=============

The ``flags`` parameter controls how pruning works:

1. Keep Only Residue
--------------------

**Flag:** ``0`` (default, no flags)

**Behavior:**

- Remove all turning points except residue
- Most aggressive pruning
- Minimal memory usage

**Use Case:** Maximum memory savings, no need for historical TPs.

**Example:**

.. code-block:: c

   // Remove all non-residue turning points
   RFC_tp_prune(&ctx, 0, 0);

   // Only residue remains
   printf("TPs after prune: %zu\n", ctx.tp_cnt);  // = residue count

2. Keep Residue + History
--------------------------

**Flag:** ``RFC_FLAGS_TPPRUNE_PRESERVE_POS``

**Behavior:**

- Keep residue
- Keep ``limit`` most recent turning points before residue
- Preserves position information

**Use Case:** Need some historical context, reprocessing capability.

**Example:**

.. code-block:: c

   // Keep residue + 100 most recent TPs
   RFC_tp_prune(&ctx, 100, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

**Result:**

.. code-block:: text

   Before: [TP1, TP2, ..., TP500, residue (50 TPs)]
   After:  [TP400, ..., TP500, residue (50 TPs)]
           ↑ 100 recent TPs     ↑ residue preserved

   Total kept: 100 + 50 = 150 TPs

How Pruning Works
=================

Algorithm Overview
------------------

1. **Identify Residue**

   Determine which turning points are part of unclosed cycles (residue).

2. **Determine Keep Set**

   Based on ``limit`` and ``flags``, decide which TPs to keep:

   - Residue (always)
   - Historical TPs (if ``RFC_FLAGS_TPPRUNE_PRESERVE_POS``)

3. **Compact Array**

   Move kept TPs to front of array, overwriting removed TPs.

4. **Update Count**

   Set ``ctx->tp_cnt`` to new count.

5. **Preserve Damage** (if RFC_DH_SUPPORT)

   Turning point damage values are preserved even during pruning.

Detailed Steps
--------------

.. code-block:: text

   Step 1: Initial State
   ─────────────────────────────────────────────────────
   TPs: [TP1, TP2, ..., TP800, residue (50 TPs)]
                                 └ Starts at TP800
   ctx->tp_cnt = 850

   Step 2: Mark Residue (Always Keep)
   ─────────────────────────────────────────────────────
   Keep: TP800 onwards (residue)

   Step 3: Add History (if PRESERVE_POS)
   ─────────────────────────────────────────────────────
   limit = 100
   Keep: TP700-TP799 (100 TPs) + residue

   Step 4: Compact Array
   ─────────────────────────────────────────────────────
   Move TP700-TP850 to beginning
   TPs: [TP700, TP701, ..., TP850, <garbage>]

   Step 5: Update Count
   ─────────────────────────────────────────────────────
   ctx->tp_cnt = 150 (100 history + 50 residue)

Automatic Pruning
=================

Instead of manually calling ``RFC_tp_prune()``, enable automatic pruning:

Configuration
-------------

.. code-block:: c

   bool RFC_tp_init_autoprune(
       void *ctx,
       bool autoprune,          // Enable/disable autopruning
       size_t size,             // Target size after prune
       size_t threshold         // Trigger threshold
   );

**Parameters:**

- ``autoprune``: ``true`` to enable, ``false`` to disable
- ``size``: Desired TP count after autopruning (passed to ``RFC_tp_prune()``)
- ``threshold``: Trigger autopruning when ``tp_cnt > threshold``

**Behavior:**

When ``ctx->tp_cnt`` exceeds ``threshold``, automatically call ``RFC_tp_prune(ctx, size, RFC_FLAGS_TPPRUNE_PRESERVE_POS)``.

Example
-------

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Enable autopruning:
   // - Trigger when TP count > 1000
   // - Prune to 500 TPs
   RFC_tp_init_autoprune(&ctx, true, 500, 1000);

   // Feed data...
   for (int i = 0; i < 1000000; i++)
   {
       RFC_feed(&ctx, &data[i], 1);

       // Autopruning happens automatically when tp_cnt > 1000
   }

**Autopruning Trigger:**

.. code-block:: text

   TP Count Over Time:
       │
   1000│        ┌┐      ┌┐         ┌┐      (Threshold)
       │   ┌────┘│  ┌───┘│  ┌──────┘│
    500│   │     └──┘    └──┘       └───── (After prune)
       ├───┘     ↑       ↑          ↑
       └─────────┴───────┴──────────┴──────> Time
                Prune    Prune    Prune

Use Cases
=========

1. Real-Time Monitoring
-----------------------

**Scenario:** Continuous monitoring system running for days/weeks.

**Solution:**

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   // Autoprun every 1000 TPs to 200
   RFC_tp_init_autoprune(&ctx, true, 200, 1000);

   // Continuous data feed
   while (monitoring_active)
   {
       double sample = read_sensor();
       RFC_feed(&ctx, &sample, 1);

       // Memory usage stays bounded!
   }

**Benefits:**

- Constant memory usage
- No manual intervention
- System runs indefinitely

2. Embedded Systems
-------------------

**Scenario:** Microcontroller with 32 KB RAM, storing TPs for offline analysis.

**Solution:**

.. code-block:: c

   // Limited RAM: keep only 100 TPs + residue
   void process_chunk(double *data, size_t len)
   {
       static rfc_ctx_s ctx;
       static bool initialized = false;

       if (!initialized)
       {
           RFC_init(&ctx, 50, 2.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
           RFC_tp_init_autoprune(&ctx, true, 50, 100);
           initialized = true;
       }

       RFC_feed(&ctx, data, len);

       // Periodically prune
       if (ctx.tp_cnt > 100)
       {
           RFC_tp_prune(&ctx, 50, RFC_FLAGS_TPPRUNE_PRESERVE_POS);
       }
   }

3. Large Dataset Processing
----------------------------

**Scenario:** Processing 100 million sample measurement file.

**Solution:**

.. code-block:: c

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);
   RFC_tp_init_autoprune(&ctx, true, 1000, 5000);

   // Process in chunks
   size_t chunk_size = 100000;
   for (size_t i = 0; i < total_samples; i += chunk_size)
   {
       size_t len = (i + chunk_size > total_samples)
           ? total_samples - i
           : chunk_size;

       RFC_feed(&ctx, &data[i], len);

       // Autopruning keeps memory bounded
       printf("TP count: %zu\n", ctx.tp_cnt);  // Never > 5000
   }

   RFC_finalize(&ctx, RFC_RES_REPEATED);

4. Circular Buffer Simulation
------------------------------

**Scenario:** Keep only last N turning points like a circular buffer.

**Solution:**

.. code-block:: c

   #define MAX_TPS 500

   void add_turning_point(rfc_ctx_s *ctx, double value)
   {
       RFC_feed(ctx, &value, 1);

       // Prune if exceeded max
       if (ctx->tp_cnt > MAX_TPS)
       {
           RFC_tp_prune(ctx, MAX_TPS, RFC_FLAGS_TPPRUNE_PRESERVE_POS);
       }
   }

Implementation Details
======================

Memory Management
-----------------

**What Happens to Removed TPs:**

.. code-block:: c

   // Before prune
   ctx->tp = [TP1, TP2, ..., TP1000]  // 1000 TPs allocated
   ctx->tp_cnt = 1000
   ctx->tp_cap = 1000

   // After prune to 100
   ctx->tp = [TP900, ..., TP1000, <garbage>]  // Same allocation
   ctx->tp_cnt = 100                            // Count reduced
   ctx->tp_cap = 1000                           // Capacity unchanged

**Note:** Pruning does NOT free memory (no realloc). It only reduces ``tp_cnt``. The allocation remains the same size.

**To Actually Free Memory:**

.. code-block:: c

   RFC_tp_prune(&ctx, 100, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Optionally shrink allocation
   ctx.tp_cap = ctx.tp_cnt;
   ctx.tp = realloc(ctx.tp, ctx.tp_cap * sizeof(rfc_value_tuple_s));

Position Field Handling
-----------------------

Each turning point has a ``pos`` field (position in original data stream):

.. code-block:: c

   struct rfc_value_tuple
   {
       double value;
       unsigned cls;
       size_t pos;      // Position in input (1-based)
       double damage;   // If RFC_DH_SUPPORT
   };

**With RFC_FLAGS_TPPRUNE_PRESERVE_POS:**

.. code-block:: c

   // Before prune
   TPs: [TP(pos=1), TP(pos=2), ..., TP(pos=1000)]

   // After prune to 100 (keep last 100)
   TPs: [TP(pos=900), TP(pos=901), ..., TP(pos=1000)]

   // Position values preserved correctly

**Without RFC_FLAGS_TPPRUNE_PRESERVE_POS:**

.. code-block:: c

   // After prune (only residue)
   TPs: [residue TPs]

   // Positions still valid, but most history lost

Damage History Interaction
---------------------------

When ``RFC_DH_SUPPORT`` is enabled, pruning interacts with damage history:

**Behavior:**

- Damage accumulated at turning points (``tp[i].damage``) is **preserved** during pruning
- Damage history array (``ctx->dh``) is **not pruned** - it tracks entire input history

**Example:**

.. code-block:: c

   // With damage history enabled
   RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, NULL, 1000000, false);

   // Feed 1,000,000 samples
   RFC_feed(&ctx, data, 1000000);

   // Prune TPs (keeps 100)
   RFC_tp_prune(&ctx, 100, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Damage history still has 1,000,000 entries
   const double *dh;
   size_t dh_count;
   RFC_dh_get(&ctx, &dh, &dh_count);
   printf("DH count: %zu\n", dh_count);  // Still 1,000,000

   // Only TP storage was pruned
   printf("TP count: %zu\n", ctx.tp_cnt);  // Now 100

Performance Considerations
==========================

Computational Cost
------------------

Pruning cost depends on how many TPs are kept:

.. code-block:: text

   Operation: Move 'keep' TPs to front of array

   Cost: O(keep_count)
         = O(limit + residue_count)

   Example:
   - Total TPs: 10,000
   - Keep: 100 + 50 residue = 150
   - Cost: O(150) ≈ 150 memory copies

   Very fast! (~microseconds)

Frequency vs. Cost Trade-off
-----------------------------

.. code-block:: text

   Strategy                     Prune Freq    Avg Memory    Total Cost
   ───────────────────────────────────────────────────────────────────
   Never prune                  Never         Maximum       0
   Prune every 10k TPs          Low           High          Low
   Prune every 1k TPs (auto)    Medium        Medium        Medium
   Prune every 100 TPs          High          Low           High

**Recommendation:** Use autopruning with threshold = 2-5× size.

.. code-block:: c

   // Good balance
   RFC_tp_init_autoprune(&ctx, true, 500, 2000);  // Prune when 2k, keep 500

Best Practices
==============

1. Enable Autopruning
---------------------

Don't call ``RFC_tp_prune()`` manually unless you have special needs:

.. code-block:: c

   // Good: Automatic
   RFC_tp_init_autoprune(&ctx, true, 500, 1000);

   // Avoid: Manual
   if (ctx.tp_cnt > 1000)
   {
       RFC_tp_prune(&ctx, 500, RFC_FLAGS_TPPRUNE_PRESERVE_POS);
   }

2. Choose Appropriate Threshold
--------------------------------

**Rule of Thumb:** ``threshold = 2 × size`` to ``5 × size``

.. code-block:: c

   // Good: 2x ratio
   RFC_tp_init_autoprune(&ctx, true, 1000, 2000);

   // Bad: Too frequent
   RFC_tp_init_autoprune(&ctx, true, 1000, 1001);  // Prunes constantly!

   // Bad: Too infrequent
   RFC_tp_init_autoprune(&ctx, true, 1000, 100000);  // Huge memory spikes

3. Balance Memory vs. Reprocessing Needs
-----------------------------------------

If you need to reprocess with different parameters later:

.. code-block:: c

   // Keep more history
   RFC_tp_init_autoprune(&ctx, true, 2000, 5000);  // Keep 2k TPs

If memory is tight and no reprocessing needed:

.. code-block:: c

   // Minimal history
   RFC_tp_init_autoprune(&ctx, true, 100, 500);  // Keep 100 TPs

4. Prune Before Finalize
-------------------------

For very long measurements, prune before calling ``RFC_finalize()``:

.. code-block:: c

   // Feed all data
   RFC_feed(&ctx, data, data_len);

   // Prune to manageable size
   RFC_tp_prune(&ctx, 1000, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Finalize (processes remaining residue)
   RFC_finalize(&ctx, RFC_RES_REPEATED);

Troubleshooting
===============

Common Issues
-------------

**1. Residue larger than limit**

.. code-block:: text

   Problem: Residue has 200 TPs, but limit = 100.

**Behavior:** Pruning keeps all residue (200 TPs), ignores limit.

**Solution:** This is correct! Residue can never be pruned. Increase limit if needed:

.. code-block:: c

   // Ensure limit > expected max residue size
   RFC_tp_prune(&ctx, 500, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

**2. TP count increases after prune**

.. code-block:: text

   Problem: tp_cnt = 100 after prune, but increases to 150 before next prune.

**Explanation:** This is normal. New turning points are added between prune calls.

.. code-block:: text

   Prune → 100 TPs → Feed data → 150 TPs → Prune → 100 TPs → ...

**3. Memory not freed**

.. code-block:: text

   Problem: Memory usage stays high after pruning.

**Cause:** Pruning doesn't call ``realloc()`` to shrink buffer.

**Solution:** Manually shrink if needed:

.. code-block:: c

   RFC_tp_prune(&ctx, 100, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Shrink allocation
   if (ctx.tp_cap > ctx.tp_cnt * 2)
   {
       ctx.tp_cap = ctx.tp_cnt;
       ctx.tp = realloc(ctx.tp, ctx.tp_cap * sizeof(rfc_value_tuple_s));
   }

**4. Position values seem wrong after prune**

.. code-block:: text

   Problem: tp[0].pos = 900, but I expected 1.

**Explanation:** Positions are preserved from original data stream. After pruning early TPs, first TP has high position value.

**This is correct!** Positions track original input, not array index.

**5. Autopruning not triggering**

.. code-block:: text

   Problem: tp_cnt grows past threshold, no pruning.

**Causes:**

- ``autoprune`` not enabled: ``RFC_tp_init_autoprune(&ctx, true, ...)``
- ``threshold`` set too high
- ``RFC_TP_SUPPORT`` not compiled in

**Verification:**

.. code-block:: c

   printf("Flags: %d\n", ctx.internal.flags & RFC_FLAGS_TPAUTOPRUNE);
   printf("Threshold: %zu\n", ctx.tp_prune_threshold);

Advanced Topics
===============

Pruning with Delegates
-----------------------

If using external TP storage (via delegates), pruning behavior changes:

.. code-block:: c

   // With external storage
   ctx.tp_set_fcn = my_tp_set;
   ctx.tp_get_fcn = my_tp_get;

   // Prune calls delegates to reorganize external storage
   RFC_tp_prune(&ctx, 100, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Your tp_set_fcn is called to rewrite TPs

See `delegates.rst <delegates.rst>`_ and `turning_points.rst <turning_points.rst>`_ for details.

Pruning and RFC_tp_refeed()
----------------------------

After pruning, you can still refeed remaining TPs:

.. code-block:: c

   // Feed initial data
   RFC_feed(&ctx, data, 10000);

   // Prune to 1000 TPs
   RFC_tp_prune(&ctx, 1000, RFC_FLAGS_TPPRUNE_PRESERVE_POS);

   // Refeed remaining TPs with different parameters
   RFC_tp_refeed(&ctx, new_hysteresis, new_class_param);

Only the remaining 1000 TPs are reprocessed.

Custom Pruning Logic
---------------------

For complete control, implement custom pruning:

.. code-block:: c

   void my_custom_prune(rfc_ctx_s *ctx)
   {
       // Your pruning algorithm
       // Example: Keep every 10th TP

       size_t write_idx = 0;
       for (size_t read_idx = 0; read_idx < ctx->tp_cnt; read_idx += 10)
       {
           if (read_idx < ctx->residue_cnt)
               continue;  // Never remove residue

           ctx->tp[write_idx++] = ctx->tp[read_idx];
       }

       ctx->tp_cnt = write_idx;
   }

Security Considerations
=======================

Integer Overflow Protection
---------------------------

The library protects against integer underflow in pruning logic:

.. code-block:: c

   // Protected calculation (from previous bug fix session)
   size_t first_to_keep = (tp_first_res_idx > limit)
       ? (tp_first_res_idx - limit)  // Safe
       : 0;                           // Avoid underflow

If ``limit`` > ``tp_first_res_idx``, could underflow without check.

**Modern versions** (after your bug fix session) include this protection.

See Also
========

- `turning_points.rst <turning_points.rst>`_ - Turning point storage and access
- `delegates.rst <delegates.rst>`_ - Custom TP storage with external backends
- `minimal_build.rst <minimal_build.rst>`_ - Memory-constrained environments
- `damage_history.rst <damage_history.rst>`_ - Interaction with damage history

References
==========

Turning point management and pruning strategies:

- **Embedded Systems** - Memory management in resource-constrained environments
- **Circular Buffers** - Common pattern for bounded storage
- **Stream Processing** - Algorithms for unbounded data streams

For implementation details:

- ``rainflow.h`` - Lines 564, 775-777 (RFC_tp_prune declaration and context fields)
- ``rainflow.c`` - Lines 894-950 (RFC_tp_prune implementation)
- ``rainflow.c`` - Lines 6701-6703 (Autopruning trigger logic)

**Bug Fix Reference:**

The integer underflow protection mentioned in "Security Considerations" was fixed during the previous code review session (before this documentation session).
