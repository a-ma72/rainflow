==============
Minimal Build
==============

Overview
========

``RFC_MINIMAL`` is a compile-time configuration that produces a minimal-footprint version of the rainflow library, specifically designed for embedded systems, microcontrollers, and resource-constrained environments. It strips away advanced features to minimize code size and memory usage while preserving core rainflow counting functionality.

What is RFC_MINIMAL?
====================

Purpose
-------

The minimal build configuration creates a lean version of the library by:

1. **Disabling advanced features** - Removes optional functionality
2. **Reducing code size** - Smaller binary footprint
3. **Minimizing RAM usage** - Lower memory requirements
4. **Simplifying API** - Smaller API surface for easier integration

This makes the library suitable for systems where resources are scarce:

- 8/16/32-bit microcontrollers
- Embedded real-time systems
- IoT devices with limited memory
- Firmware with strict size constraints

What Gets Disabled?
===================

When ``RFC_MINIMAL`` is enabled, the following features are automatically disabled:

Disabled Features
-----------------

.. code-block:: c

   #if RFC_MINIMAL
   #undef  RFC_USE_DELEGATES
   #define RFC_USE_DELEGATES    OFF    // No delegate functions

   #undef  RFC_HCM_SUPPORT
   #define RFC_HCM_SUPPORT      OFF    // No HCM counting method

   #undef  RFC_ASTM_SUPPORT
   #define RFC_ASTM_SUPPORT     OFF    // No ASTM stack-based algorithm

   #undef  RFC_TP_SUPPORT
   #define RFC_TP_SUPPORT       OFF    // No turning point storage

   #undef  RFC_DH_SUPPORT
   #define RFC_DH_SUPPORT       OFF    // No damage history

   #undef  RFC_AT_SUPPORT
   #define RFC_AT_SUPPORT       OFF    // No amplitude transformation

   #undef  RFC_AR_SUPPORT
   #define RFC_AR_SUPPORT       OFF    // No amplitude-range support

   #undef  RFC_GLOBAL_EXTREMA
   #define RFC_GLOBAL_EXTREMA   OFF    // No global min/max tracking

   #undef  RFC_DAMAGE_FAST
   #define RFC_DAMAGE_FAST      OFF    // No damage lookup tables
   #endif

**Summary Table:**

+-------------------------+------------------+-----------------------------+
| Feature                 | Status           | Impact                      |
+=========================+==================+=============================+
| RFC_USE_DELEGATES       | Disabled         | No custom callbacks         |
+-------------------------+------------------+-----------------------------+
| RFC_HCM_SUPPORT         | Disabled         | Only 4-point method         |
+-------------------------+------------------+-----------------------------+
| RFC_ASTM_SUPPORT        | Disabled         | No ASTM stack algorithm     |
+-------------------------+------------------+-----------------------------+
| RFC_TP_SUPPORT          | Disabled         | No TP storage/access        |
+-------------------------+------------------+-----------------------------+
| RFC_DH_SUPPORT          | Disabled         | No damage history           |
+-------------------------+------------------+-----------------------------+
| RFC_AT_SUPPORT          | Disabled         | No Haigh diagrams           |
+-------------------------+------------------+-----------------------------+
| RFC_AR_SUPPORT          | Disabled         | No amplitude-range pairs    |
+-------------------------+------------------+-----------------------------+
| RFC_GLOBAL_EXTREMA      | Disabled         | No global min/max           |
+-------------------------+------------------+-----------------------------+
| RFC_DAMAGE_FAST         | Disabled         | No damage LUT optimization  |
+-------------------------+------------------+-----------------------------+

What Remains Enabled?
======================

Core Functionality
------------------

The minimal build retains all essential rainflow counting features:

**Available Features:**

1. **Basic Rainflow Counting**

   - 4-point rainflow algorithm
   - Cycle detection and counting
   - Rainflow matrix generation

2. **Damage Calculation**

   - Wöhler curve (S-N curve) support
   - Miner's rule damage accumulation
   - Elementary and modified Wöhler curves

3. **Histograms**

   - Rainflow matrix (from-to)
   - Range-pair counting
   - Level crossing counting (basic)

4. **Residue Handling** (Limited)

   - ``RFC_RES_NONE`` - Ignore residue
   - ``RFC_RES_IGNORE`` - Same as NONE
   - ``RFC_RES_NO_FINALIZE`` - Don't finalize

   **Not Available:**

   - RFC_RES_HALFCYCLES
   - RFC_RES_FULLCYCLES
   - RFC_RES_CLORMANN_SEEGER
   - RFC_RES_REPEATED
   - RFC_RES_RP_DIN45667

5. **Basic Configuration**

   - Class width/offset/count
   - Hysteresis filtering (if enabled separately)

Memory Footprint Comparison
============================

Typical Size Comparison
-----------------------

Approximate code sizes (compiled for ARM Cortex-M4, -Os optimization):

.. code-block:: text

   Configuration          Code Size    Data/BSS    Total ROM
   ──────────────────────────────────────────────────────────
   Full Build (all)       ~45 KB       ~2 KB       ~47 KB
   Typical Build          ~35 KB       ~1.5 KB     ~36.5 KB
   RFC_MINIMAL            ~15 KB       ~0.5 KB     ~15.5 KB

   Reduction: ~65% code size savings

Runtime RAM Usage
-----------------

For 100 classes, typical data:

.. code-block:: text

   Configuration          rfc_ctx_s    Histograms  Total RAM
   ──────────────────────────────────────────────────────────
   Full Build             ~500 bytes   80 KB       ~81 KB
   RFC_MINIMAL            ~200 bytes   80 KB       ~80 KB

   (Histogram size dominates, context size reduced by ~60%)

Enabling RFC_MINIMAL
====================

Compile-Time Configuration
---------------------------

**Method 1: CMake**

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_MINIMAL=ON

**Method 2: config.h**

.. code-block:: c

   // In config.h or build system
   #define RFC_MINIMAL 1

**Method 3: Preprocessor Define**

.. code-block:: bash

   gcc -DRFC_MINIMAL=1 -c rainflow.c

**Verification:**

Check that RFC_MINIMAL is active:

.. code-block:: c

   #include "rainflow.h"

   #if RFC_MINIMAL
   printf("RFC_MINIMAL is enabled\n");
   #else
   printf("RFC_MINIMAL is disabled\n");
   #endif

Minimal Build API
=================

Available Functions
-------------------

The following API functions are available in RFC_MINIMAL mode:

**Initialization:**

.. code-block:: c

   bool RFC_init(rfc_ctx_s *ctx,
                 unsigned class_count,
                 double class_width,
                 double class_offset,
                 double hysteresis,
                 rfc_flags_e flags);

   bool RFC_deinit(rfc_ctx_s *ctx);

**Data Processing:**

.. code-block:: c

   bool RFC_feed(rfc_ctx_s *ctx, const double *data, size_t data_len);
   bool RFC_finalize(rfc_ctx_s *ctx, rfc_res_method_e residual_method);

**Wöhler Curve:**

.. code-block:: c

   bool RFC_wl_init_elementary(rfc_ctx_s *ctx, double sx, double nx, double k);
   bool RFC_wl_init_modified(rfc_ctx_s *ctx, double sx, double nx, double k, double k2);

**Results:**

.. code-block:: c

   bool RFC_damage(rfc_ctx_s *ctx, double *damage);
   bool RFC_rfm_get(rfc_ctx_s *ctx, const double **rfm, size_t *count);
   bool RFC_rp_get(rfc_ctx_s *ctx, const double **rp, size_t *count);
   bool RFC_lc_get(rfc_ctx_s *ctx, const double **lc, size_t *count);
   bool RFC_res_get(rfc_ctx_s *ctx, const rfc_value_tuple_s **residue, unsigned *count);

**NOT Available:**

- RFC_tp_* functions (turning point access)
- RFC_dh_* functions (damage history)
- RFC_at_* functions (amplitude transformation)
- Delegate-related functions

Example Usage
=============

Complete Minimal Example
-------------------------

.. code-block:: c

   #include "rainflow.h"
   #include <stdio.h>

   int main(void)
   {
       rfc_ctx_s ctx;
       double data[] = {0, 10, 0, 20, 0, 30, 0, 20, 10, 0};
       size_t data_len = sizeof(data) / sizeof(data[0]);

       // Initialize with 100 classes, 1.0 class width
       if (!RFC_init(&ctx, 100, 1.0, 0.0, 0.5, RFC_FLAGS_DEFAULT))
       {
           fprintf(stderr, "RFC_init failed\n");
           return 1;
       }

       // Set up Wöhler curve for damage calculation
       RFC_wl_init_elementary(&ctx, 1000.0, 1e7, 5.0);

       // Feed data
       if (!RFC_feed(&ctx, data, data_len))
       {
           fprintf(stderr, "RFC_feed failed\n");
           RFC_deinit(&ctx);
           return 1;
       }

       // Finalize (RFC_MINIMAL only supports NONE/IGNORE/NO_FINALIZE)
       if (!RFC_finalize(&ctx, RFC_RES_IGNORE))
       {
           fprintf(stderr, "RFC_finalize failed\n");
           RFC_deinit(&ctx);
           return 1;
       }

       // Get results
       double damage;
       RFC_damage(&ctx, &damage);
       printf("Total damage: %.6e\n", damage);

       printf("Full cycles: %.1f\n", ctx.full_inc);
       printf("Half cycles: %.1f\n", ctx.half_inc);

       // Get rainflow matrix
       const double *rfm;
       size_t rfm_count;
       RFC_rfm_get(&ctx, &rfm, &rfm_count);

       printf("\nRainflow Matrix (non-zero entries):\n");
       for (size_t i = 0; i < ctx.class_count; i++)
       {
           for (size_t j = 0; j < ctx.class_count; j++)
           {
               double count = rfm[i * ctx.class_count + j];
               if (count > 0.0)
               {
                   printf("  [%zu → %zu]: %.2f cycles\n", i, j, count);
               }
           }
       }

       // Cleanup
       RFC_deinit(&ctx);

       return 0;
   }

Embedded System Example
------------------------

For a microcontroller with severe RAM constraints:

.. code-block:: c

   // Static allocation example (no malloc)
   #include "rainflow.h"

   // Pre-allocated static buffers
   static rfc_ctx_s ctx;
   static double rfm_buffer[50 * 50];  // 50 classes = 10KB
   static rfc_value_tuple_s residue_buffer[100];

   void embedded_rainflow_init(void)
   {
       // Initialize with static buffers
       RFC_init(&ctx, 50, 2.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

       // Attach static buffers (manual setup)
       ctx.rfm = rfm_buffer;
       ctx.residue = residue_buffer;
       ctx.residue_cap = 100;

       // Wöhler curve
       RFC_wl_init_elementary(&ctx, 800.0, 1e6, 4.0);
   }

   void embedded_rainflow_process(const double *data, size_t len)
   {
       RFC_feed(&ctx, data, len);
   }

   double embedded_rainflow_finalize(void)
   {
       RFC_finalize(&ctx, RFC_RES_IGNORE);

       double damage;
       RFC_damage(&ctx, &damage);
       return damage;
   }

   // No need to call RFC_deinit() if using static allocation

Porting to Microcontrollers
============================

Target Platforms
----------------

RFC_MINIMAL is suitable for:

- **ARM Cortex-M0/M0+**: 32 KB Flash, 8 KB RAM
- **ARM Cortex-M3/M4**: 128 KB Flash, 32 KB RAM
- **AVR ATmega**: 32-256 KB Flash, 2-8 KB RAM
- **ESP32**: Plenty of resources, but minimal build reduces overhead
- **STM32**: Most STM32F0/F1/F4 series
- **PIC32**: Microchip 32-bit microcontrollers

Build Considerations
--------------------

**1. No Dynamic Allocation (Optional)**

If your system doesn't have malloc():

.. code-block:: c

   // Define custom allocators
   #define RFC_ALLOC(size) my_static_alloc(size)
   #define RFC_FREE(ptr) my_static_free(ptr)

Or use fully static allocation (see embedded example above).

**2. Fixed-Point Math (Optional)**

For systems without FPU, consider using fixed-point:

.. code-block:: c

   // Use integer types for values
   #define RFC_VALUE_TYPE int32_t  // Fixed-point Q16.16

   // Scale input data appropriately
   int32_t data_fixed[] = {0, 10 << 16, 0, 20 << 16, ...};

**3. Reduce Class Count**

Fewer classes = less memory:

.. code-block:: c

   // Instead of 100 classes (80 KB matrix)
   RFC_init(&ctx, 100, 1.0, 0.0, 0.5, RFC_FLAGS_DEFAULT);

   // Use 20 classes (3.2 KB matrix)
   RFC_init(&ctx, 20, 5.0, 0.0, 2.5, RFC_FLAGS_DEFAULT);

Memory Usage Formula
--------------------

.. code-block:: text

   RAM = ctx_size + rfm_size + residue_size

   Where:
   ctx_size     ≈ 200 bytes (RFC_MINIMAL)
   rfm_size     = class_count² × sizeof(double)
                = class_count² × 8 bytes
   residue_size = residue_cap × sizeof(rfc_value_tuple_s)
                = residue_cap × 24 bytes

   Example (50 classes, 50 residue):
   RAM = 200 + (50² × 8) + (50 × 24)
       = 200 + 20,000 + 1,200
       = 21,400 bytes ≈ 21 KB

Performance Considerations
==========================

Speed Comparison
----------------

RFC_MINIMAL is typically **10-20% faster** than full build due to:

- Fewer conditional checks
- Smaller code fits better in instruction cache
- No overhead from disabled features

Benchmark (1M samples, ARM Cortex-M4 @ 168 MHz):

.. code-block:: text

   Configuration          Processing Time    Speedup
   ────────────────────────────────────────────────
   Full Build             850 ms             1.0x
   RFC_MINIMAL            720 ms             1.18x

Optimization Tips
-----------------

**1. Compiler Optimization**

.. code-block:: bash

   # Size optimization
   gcc -Os -DRFC_MINIMAL=1 -c rainflow.c

   # Speed optimization (if Flash permits)
   gcc -O2 -DRFC_MINIMAL=1 -c rainflow.c

**2. Link-Time Optimization (LTO)**

.. code-block:: bash

   gcc -flto -Os -DRFC_MINIMAL=1 rainflow.c main.c -o firmware.elf

**3. Profile-Guided Optimization (PGO)**

For embedded systems with profiling capability, PGO can yield additional 5-10% speedup.

Limitations of RFC_MINIMAL
===========================

What You Cannot Do
-------------------

1. **No Turning Point Storage**

   Cannot access turning point history:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // RFC_tp_get(&ctx, 5, &tp);  // Compile error!

2. **No Damage History**

   Cannot track damage over time:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, ...);  // Compile error!

3. **No ASTM Method**

   Only 4-point algorithm available:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // ctx.counting_method = RFC_COUNTING_METHOD_ASTM;  // Won't work

4. **No HCM Method**

   Clormann/Seeger method unavailable:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // ctx.counting_method = RFC_COUNTING_METHOD_HCM;  // Won't work

5. **Limited Residue Methods**

   Only basic residue handling:

   .. code-block:: c

      // AVAILABLE:
      RFC_finalize(&ctx, RFC_RES_NONE);
      RFC_finalize(&ctx, RFC_RES_IGNORE);
      RFC_finalize(&ctx, RFC_RES_NO_FINALIZE);

      // NOT AVAILABLE in RFC_MINIMAL:
      // RFC_finalize(&ctx, RFC_RES_REPEATED);  // Compile error!

6. **No Delegates**

   Cannot customize behavior via function pointers:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // ctx.cycle_find_fcn = my_cycle_finder;  // Field doesn't exist

7. **No Amplitude Transformation**

   Haigh diagrams and mean stress correction unavailable:

   .. code-block:: c

      // NOT AVAILABLE in RFC_MINIMAL
      // RFC_at_init(&ctx, ...);  // Compile error!

Workarounds
-----------

**If you need disabled features:**

1. **Use selective compilation** - Enable only needed features instead of full build:

   .. code-block:: bash

      cmake -S. -Bbuild \
          -DRFC_TP_SUPPORT=OFF \
          -DRFC_DH_SUPPORT=OFF \
          -DRFC_HCM_SUPPORT=OFF \
          -DRFC_ASTM_SUPPORT=OFF

   This saves memory while keeping some features.

2. **Post-process on PC** - Collect raw data on embedded device, perform full analysis on PC.

3. **External storage** - Even without RFC_TP_SUPPORT, you can log data separately.

Migration Guide
===============

From Full Build to RFC_MINIMAL
-------------------------------

**Step 1: Identify Dependencies**

Check your code for disabled features:

.. code-block:: bash

   # Search for features not in RFC_MINIMAL
   grep -r "RFC_tp_" your_code/
   grep -r "RFC_dh_" your_code/
   grep -r "RFC_at_" your_code/
   grep -r "\.counting_method.*ASTM" your_code/

**Step 2: Remove or Conditionalize**

Wrap feature usage in preprocessor guards:

.. code-block:: c

   #if !RFC_MINIMAL
   RFC_dh_init(&ctx, RFC_SD_TRANSIENT_23c, NULL, data_len, false);
   #endif

   // Or provide alternative
   #if RFC_TP_SUPPORT
   RFC_tp_get(&ctx, 5, &tp);
   #else
   // Alternative implementation without TP storage
   #endif

**Step 3: Simplify Residue Handling**

Replace advanced residue methods:

.. code-block:: c

   // Before (full build)
   RFC_finalize(&ctx, RFC_RES_REPEATED);

   // After (RFC_MINIMAL compatible)
   RFC_finalize(&ctx, RFC_RES_IGNORE);

**Step 4: Rebuild and Test**

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_MINIMAL=ON
   cmake --build build
   # Test thoroughly!

Testing RFC_MINIMAL
===================

Validation
----------

Ensure RFC_MINIMAL produces correct results:

.. code-block:: c

   // Test against known dataset
   double test_data[] = {0, 10, 0, 20, 0, 30, 0};

   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 0.5, RFC_FLAGS_DEFAULT);
   RFC_wl_init_elementary(&ctx, 1000.0, 1e7, 5.0);

   RFC_feed(&ctx, test_data, 7);
   RFC_finalize(&ctx, RFC_RES_IGNORE);

   double damage;
   RFC_damage(&ctx, &damage);

   // Verify expected damage (compare with full build)
   const double expected_damage = 3.2e-5;  // From full build
   assert(fabs(damage - expected_damage) < 1e-10);

   RFC_deinit(&ctx);

Memory Profiling
----------------

Measure actual RAM usage on target:

.. code-block:: c

   extern unsigned int _heap_start;
   extern unsigned int _heap_end;

   void* heap_start = &_heap_start;
   void* heap_before = sbrk(0);

   rfc_ctx_s ctx;
   RFC_init(&ctx, 50, 1.0, 0.0, 0.5, RFC_FLAGS_DEFAULT);

   void* heap_after = sbrk(0);
   size_t memory_used = (size_t)(heap_after - heap_before);

   printf("Memory allocated: %zu bytes\n", memory_used);

   RFC_deinit(&ctx);

See Also
========

- `features.rst <features.rst>`_ - All compile-time feature flags
- `installation.rst <installation.rst>`_ - Build system configuration
- `delegates.rst <delegates.rst>`_ - Why delegates are disabled in RFC_MINIMAL
- `algorithm.rst <algorithm.rst>`_ - Core algorithm (always available)

References
==========

Embedded systems programming resources:

- **ARM Cortex-M Programming** - Joseph Yiu, "The Definitive Guide"
- **Embedded C** - Michael Barr, "Programming Embedded Systems in C and C++"
- **Microcontroller Optimization** - Compiler and linker optimization guides

For RFC_MINIMAL implementation details:

- ``rainflow.h`` - Lines 183-201 (RFC_MINIMAL feature disabling)
- ``rainflow.c`` - Conditional compilation blocks
- ``CMakeLists.txt`` - RFC_MINIMAL build configuration
