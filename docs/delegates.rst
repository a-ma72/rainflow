=========
Delegates
=========

Overview
========

Delegates are a powerful extensibility mechanism in the rainflow library that allow you to customize core behavior by providing your own callback functions. Think of delegates as "hooks" where you can inject custom logic into the rainflow counting process.

What are Delegates?
===================

Concept
-------

**Delegates** are function pointers that the rainflow library calls at specific points during processing. When set, your custom function is called instead of (or in addition to) the default implementation.

This allows you to:

- **Customize counting algorithms** - Implement your own cycle detection logic
- **Control memory management** - Use custom allocators or external storage
- **Modify damage calculations** - Apply custom fatigue models
- **Transform inputs** - Apply custom amplitude transformations
- **Debug and log** - Intercept internal operations for monitoring

Key Benefits
------------

**1. Flexibility Without Forking**

   Customize behavior without modifying library source code.

**2. Memory Control**

   Integrate with custom memory management systems (pools, arenas, etc.).

**3. External Storage**

   Store turning points in databases, files, or specialized data structures.

**4. Custom Algorithms**

   Implement proprietary counting methods or material models.

**5. Real-Time Integration**

   Hook into embedded systems with specific constraints.

Enabling Delegates
==================

Compile-Time Configuration
---------------------------

Delegates must be enabled during compilation:

.. code-block:: bash

   cmake -S. -Bbuild -DRFC_USE_DELEGATES=ON

Or in your config:

.. code-block:: c

   #define RFC_USE_DELEGATES 1

**Important:** RFC_MINIMAL automatically disables delegates to minimize code size.

Available Delegates
===================

The library provides several delegate function types. All delegates are optional - if not set (NULL), the library uses its default implementation.

1. Cycle Finding Delegate
--------------------------

**Type:** ``rfc_cycle_find_fcn_t``

**Signature:**

.. code-block:: c

   typedef void (*rfc_cycle_find_fcn_t)(rfc_ctx_s *rfc_ctx, rfc_flags_e flags);

**Purpose:** Replace the entire cycle counting algorithm with your own implementation.

**When Called:** Every time ``RFC_cycle_find()`` is invoked during rainflow counting.

**Use Cases:**

- Implement proprietary counting algorithms
- Integrate third-party cycle detection methods
- Apply custom filtering or preprocessing

**Example:**

.. code-block:: c

   void my_cycle_finder(rfc_ctx_s *ctx, rfc_flags_e flags)
   {
       // Custom cycle detection logic
       // Access turning points via ctx->tp
       // Count cycles by calling RFC_count_cycle()

       printf("Custom cycle finder called with %zu turning points\n",
              ctx->tp_cnt);

       // Your algorithm here...
   }

   // Usage
   rfc_ctx_s ctx;
   RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

   ctx.cycle_find_fcn = my_cycle_finder;
   ctx.counting_method = RFC_COUNTING_METHOD_DELEGATED;  // Important!

2. Damage Calculation Delegate
-------------------------------

**Type:** ``rfc_damage_calc_fcn_t``

**Signature:**

.. code-block:: c

   typedef bool (*rfc_damage_calc_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       unsigned from_class,    // Starting class
       unsigned to_class,      // Ending class
       double *damage,         // [out] Computed damage
       double *Sa_ret          // [out] Stress amplitude
   );

**Purpose:** Override damage calculation for each cycle.

**When Called:** Every time a cycle is counted and damage needs to be calculated.

**Use Cases:**

- Apply custom Wöhler curves
- Implement non-standard fatigue models
- Integrate material database lookups
- Apply temperature-dependent damage

**Example:**

.. code-block:: c

   bool my_damage_calculator(rfc_ctx_s *ctx,
                             unsigned from_class,
                             unsigned to_class,
                             double *damage,
                             double *Sa_ret)
   {
       // Calculate stress amplitude
       double range = fabs((double)to_class - (double)from_class) * ctx->class_width;
       double Sa = range / 2.0;

       // Custom damage model (e.g., Coffin-Manson for strain-based)
       double epsilon_a = Sa / YOUNGS_MODULUS;  // Your material property
       double Nf = pow(epsilon_a / 0.01, -0.6); // Your model

       *damage = 1.0 / Nf;
       *Sa_ret = Sa;

       return true;  // Success
   }

   // Usage
   ctx.damage_calc_fcn = my_damage_calculator;

3. Finalize Delegate
--------------------

**Type:** ``rfc_finalize_fcn_t``

**Signature:**

.. code-block:: c

   typedef bool (*rfc_finalize_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       int residual_method     // RFC_RES_* method
   );

**Purpose:** Custom residue handling logic.

**When Called:** When ``RFC_finalize()`` is called to process residue.

**Use Cases:**

- Implement custom residue methods
- Apply special handling for unclosed cycles
- Skip residue processing entirely

**Example:**

.. code-block:: c

   bool my_finalize(rfc_ctx_s *ctx, int residual_method)
   {
       printf("Custom finalize: %zu residue points\n", ctx->residue_cnt);

       // Custom residue handling
       if (residual_method == RFC_RES_IGNORE)
       {
           // Do nothing with residue
           return true;
       }

       // Custom processing...
       return true;
   }

   // Usage
   ctx.finalize_fcn = my_finalize;

4. Turning Point Delegates (with RFC_TP_SUPPORT)
-------------------------------------------------

When ``RFC_TP_SUPPORT`` is enabled, three delegates control turning point storage:

**A. tp_next_fcn - Turning Point Detection**

**Type:** ``rfc_tp_next_fcn_t``

**Signature:**

.. code-block:: c

   typedef rfc_value_tuple_s* (*rfc_tp_next_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       const rfc_value_tuple_s *pt
   );

**Purpose:** Decide if a new value forms a turning point.

**When Called:** For each new input value during ``RFC_feed()``.

**Example:**

.. code-block:: c

   rfc_value_tuple_s* my_tp_next(rfc_ctx_s *ctx, const rfc_value_tuple_s *pt)
   {
       // Apply custom turning point detection
       // Return pt if it's a turning point, NULL otherwise

       // Example: Add hysteresis filtering
       double min_change = 10.0;  // Minimum significant change

       if (fabs(pt->value - ctx->internal.pos.value) > min_change)
       {
           return (rfc_value_tuple_s*)pt;  // Accept as turning point
       }

       return NULL;  // Reject - not significant enough
   }

   ctx.tp_next_fcn = my_tp_next;

**B. tp_set_fcn - Write Turning Point**

**Type:** ``rfc_tp_set_fcn_t``

**Signature:**

.. code-block:: c

   typedef bool (*rfc_tp_set_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       size_t tp_pos,           // Position (1-based)
       rfc_value_tuple_s *tp    // Turning point data
   );

**Purpose:** Store a turning point in custom storage.

**When Called:** When a turning point needs to be written.

**Use Cases:**

- Store turning points in a database
- Write to external memory or file
- Implement circular buffer storage
- Log to network endpoint

**Example:**

.. code-block:: c

   bool my_tp_set(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
   {
       // Write to custom storage
       // tp_pos is 1-based

       // Example: Write to file
       FILE *fp = (FILE*)ctx->internal.obj;  // Custom storage handle
       if (fp)
       {
           fprintf(fp, "%zu,%f,%u,%zu\n",
                   tp_pos, tp->value, tp->cls, tp->pos);
           return true;
       }

       return false;
   }

   // Usage
   FILE *tp_file = fopen("turning_points.csv", "w");
   ctx.internal.obj = tp_file;  // Store handle
   ctx.tp_set_fcn = my_tp_set;

**C. tp_get_fcn - Read Turning Point**

**Type:** ``rfc_tp_get_fcn_t``

**Signature:**

.. code-block:: c

   typedef bool (*rfc_tp_get_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       size_t tp_pos,               // Position (1-based)
       rfc_value_tuple_s **tp       // [out] Pointer to turning point
   );

**Purpose:** Retrieve a turning point from custom storage.

**When Called:** When the algorithm needs to read a turning point.

**Example:**

.. code-block:: c

   // External storage
   rfc_value_tuple_s *external_tp_storage = NULL;
   size_t external_tp_capacity = 0;

   bool my_tp_get(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s **tp)
   {
       if (!tp_pos || tp_pos > external_tp_capacity)
           return false;

       // Return pointer to turning point in external storage
       *tp = &external_tp_storage[tp_pos - 1];  // tp_pos is 1-based
       return true;
   }

   // Setup
   external_tp_capacity = 10000;
   external_tp_storage = malloc(external_tp_capacity * sizeof(rfc_value_tuple_s));
   ctx.tp_get_fcn = my_tp_get;
   ctx.tp_set_fcn = my_tp_set;  // Need both get and set

5. Damage History Delegate (with RFC_DH_SUPPORT)
-------------------------------------------------

**Type:** ``rfc_spread_damage_fcn_t``

**Signature:**

.. code-block:: c

   typedef void (*rfc_spread_damage_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       rfc_value_tuple_s *from,     // Cycle start
       rfc_value_tuple_s *to,       // Cycle end
       rfc_value_tuple_s *next,     // Next point
       rfc_flags_e flags
   );

**Purpose:** Custom damage distribution logic for damage history.

**When Called:** When a cycle closes and damage needs to be spread across the timeline.

**Use Cases:**

- Implement custom transient damage models
- Apply physics-based damage spreading
- Integrate with FEM solvers

See `spread_damage.rst <spread_damage.rst>`_ for detailed explanation.

6. Amplitude Transformation Delegate (with RFC_AT_SUPPORT)
-----------------------------------------------------------

**Type:** ``rfc_at_transform_fcn_t``

**Signature:**

.. code-block:: c

   typedef bool (*rfc_at_transform_fcn_t)(
       rfc_ctx_s *rfc_ctx,
       double Sa,                    // Stress amplitude
       double Sm,                    // Mean stress
       double *Sa_transformed        // [out] Transformed amplitude
   );

**Purpose:** Apply mean stress correction (Haigh diagram) with custom model.

**When Called:** During damage calculation when amplitude transformation is active.

**Use Cases:**

- Implement custom mean stress corrections
- Apply material-specific transformations
- Use proprietary damage models

**Example:**

.. code-block:: c

   bool my_haigh_transform(rfc_ctx_s *ctx, double Sa, double Sm,
                           double *Sa_transformed)
   {
       // Custom mean stress correction
       // Example: Modified Goodman equation
       double Su = 1000.0;  // Ultimate tensile strength
       double correction = 1.0 - (Sm / Su);

       if (correction <= 0.0)
       {
           *Sa_transformed = INFINITY;  // Failure region
       }
       else
       {
           *Sa_transformed = Sa / correction;
       }

       return true;
   }

   ctx.at_transform_fcn = my_haigh_transform;

7. Debug Delegate (with RFC_DEBUG_FLAGS)
-----------------------------------------

**Type:** ``rfc_debug_vfprintf_fcn_t``

**Signature:**

.. code-block:: c

   typedef int (*rfc_debug_vfprintf_fcn_t)(
       void *obj,
       FILE *stream,
       const char *fmt,
       va_list arg
   );

**Purpose:** Intercept debug output for custom logging.

**When Called:** When library outputs debug information.

**Use Cases:**

- Redirect debug output to custom log system
- Filter or format debug messages
- Send logs to network or database

**Example:**

.. code-block:: c

   int my_debug_logger(void *obj, FILE *stream, const char *fmt, va_list arg)
   {
       // Custom logging
       char buffer[1024];
       int result = vsnprintf(buffer, sizeof(buffer), fmt, arg);

       // Send to custom log system
       my_log_system_write("RFC", buffer);

       return result;
   }

   ctx.debug_vfprintf_fcn = my_debug_logger;

Complete Example: Custom Storage Backend
=========================================

Here's a complete example showing how to use delegates for external turning point storage:

.. code-block:: c

   #include "rainflow.h"
   #include <stdio.h>
   #include <stdlib.h>

   // External storage structure
   typedef struct {
       FILE *file;
       rfc_value_tuple_s *cache;
       size_t cache_capacity;
   } tp_storage_s;

   // Set turning point (write to file)
   bool tp_set_file(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
   {
       tp_storage_s *storage = (tp_storage_s*)ctx->internal.obj;

       // Also keep in memory cache
       if (tp_pos > 0 && tp_pos <= storage->cache_capacity)
       {
           storage->cache[tp_pos - 1] = *tp;
       }

       // Write to file
       if (storage->file)
       {
           fseek(storage->file, (tp_pos - 1) * sizeof(rfc_value_tuple_s), SEEK_SET);
           fwrite(tp, sizeof(rfc_value_tuple_s), 1, storage->file);
           fflush(storage->file);
       }

       return true;
   }

   // Get turning point (read from cache or file)
   bool tp_get_file(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s **tp)
   {
       tp_storage_s *storage = (tp_storage_s*)ctx->internal.obj;

       if (!tp_pos || tp_pos > storage->cache_capacity)
           return false;

       // Return from cache
       *tp = &storage->cache[tp_pos - 1];

       return true;
   }

   int main(void)
   {
       rfc_ctx_s ctx;
       tp_storage_s storage;
       double data[] = {0, 10, 0, 20, 0, 30, 0};

       // Initialize storage
       storage.file = fopen("turning_points.bin", "wb+");
       storage.cache_capacity = 1000;
       storage.cache = calloc(storage.cache_capacity, sizeof(rfc_value_tuple_s));

       // Initialize rainflow with delegates
       RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT);

       // Attach storage to context
       ctx.internal.obj = &storage;

       // Set delegates
       ctx.tp_set_fcn = tp_set_file;
       ctx.tp_get_fcn = tp_get_file;

       // Process data
       RFC_feed(&ctx, data, sizeof(data) / sizeof(data[0]));
       RFC_finalize(&ctx, RFC_RES_REPEATED);

       // Print results
       double damage;
       RFC_damage(&ctx, &damage);
       printf("Total damage: %.6e\n", damage);

       // Cleanup
       RFC_deinit(&ctx);
       fclose(storage.file);
       free(storage.cache);

       return 0;
   }

Using the internal.obj Field
=============================

The ``rfc_ctx_s`` structure includes a ``void *obj`` field specifically for delegate use:

.. code-block:: c

   struct rfc_ctx_s
   {
       // ... other fields ...

       #if RFC_USE_DELEGATES
       struct {
           void *obj;  // User-defined object for delegates
       } internal;
       #endif
   };

**Purpose:** Store custom context data that your delegates can access.

**Usage:**

.. code-block:: c

   // Store custom data
   my_custom_data_s my_data = { /* ... */ };
   ctx.internal.obj = &my_data;

   // Access in delegate
   bool my_delegate(rfc_ctx_s *ctx, /* ... */)
   {
       my_custom_data_s *data = (my_custom_data_s*)ctx->internal.obj;
       // Use data...
   }

Delegate Best Practices
=======================

1. Check for NULL
-----------------

Always check if a delegate is set before calling:

.. code-block:: c

   #if RFC_USE_DELEGATES
   if (ctx->tp_set_fcn)
   {
       ctx->tp_set_fcn(ctx, tp_pos, tp);
   }
   else
   #endif
   {
       // Default implementation
   }

The library does this internally.

2. Return Success/Failure
--------------------------

Delegates that return ``bool`` should:

- Return ``true`` on success
- Return ``false`` on error

The library will handle errors appropriately.

3. Thread Safety
----------------

Delegates are called from the same thread as ``RFC_feed()`` and ``RFC_finalize()``. If your delegate accesses shared resources, ensure proper synchronization.

4. Memory Management
--------------------

If your delegate allocates memory, ensure it's freed:

- Store allocations in ``ctx->internal.obj``
- Free during ``RFC_deinit()`` or in a cleanup function

5. Performance Considerations
------------------------------

Delegates are called frequently (potentially millions of times):

- Keep delegate functions fast
- Avoid I/O in performance-critical delegates
- Use buffering for logging/storage operations

C++ Integration
===============

Delegates work seamlessly with C++ by using lambda captures or member functions:

Using Lambdas
-------------

.. code-block:: cpp

   #include "rainflow.hpp"

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);

   // Use lambda as delegate
   auto my_damage_calc = [](rfc_ctx_s *ctx, unsigned from_class,
                           unsigned to_class, double *damage, double *Sa)
   {
       // Custom calculation
       *damage = 0.001;
       *Sa = 100.0;
       return true;
   };

   rf.ctx_get().damage_calc_fcn = my_damage_calc;

Using Member Functions
----------------------

.. code-block:: cpp

   class CustomRainflow
   {
   public:
       Rainflow::Rainflow rf;

       CustomRainflow()
       {
           rf.init(100, 1.0);

           // Store 'this' pointer
           rf.ctx_get().internal.obj = this;

           // Set delegate
           rf.ctx_get().damage_calc_fcn = &CustomRainflow::damage_calc_static;
       }

       // Static wrapper
       static bool damage_calc_static(rfc_ctx_s *ctx, unsigned from_class,
                                     unsigned to_class, double *damage, double *Sa)
       {
           CustomRainflow *self = (CustomRainflow*)ctx->internal.obj;
           return self->damage_calc_member(from_class, to_class, damage, Sa);
       }

       // Member function with access to class members
       bool damage_calc_member(unsigned from_class, unsigned to_class,
                              double *damage, double *Sa)
       {
           // Can access this->material_properties, etc.
           *damage = calculate_damage_from_material_db(from_class, to_class);
           return true;
       }

   private:
       // Material properties, databases, etc.
   };

Performance Impact
==================

Overhead
--------

When delegates are:

- **Not set (NULL)**: Zero overhead - delegates are checked with simple pointer comparison
- **Set**: One function pointer call overhead (~1-5 CPU cycles)

For most applications, delegate overhead is negligible compared to the computational cost of counting and damage calculation.

Optimization Tips
-----------------

**1. Inline Small Delegates**

   Delegates are called through function pointers, so they won't be inlined. Keep them simple.

**2. Batch Operations**

   If your delegate does I/O, buffer operations:

   .. code-block:: c

      // Bad: Write each turning point immediately
      bool tp_set_slow(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
      {
          fwrite(tp, sizeof(*tp), 1, file);  // Slow!
          return true;
      }

      // Good: Buffer writes
      bool tp_set_fast(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
      {
          buffer[buffer_count++] = *tp;
          if (buffer_count >= BUFFER_SIZE)
          {
              fwrite(buffer, sizeof(*tp), buffer_count, file);
              buffer_count = 0;
          }
          // TODO: Flush remaining buffer on cleanup
          return true;
      }

**3. Avoid Dynamic Allocation**

   Pre-allocate resources during setup, not in delegates.

Use Cases and Examples
======================

1. Embedded Systems with Limited RAM
-------------------------------------

**Problem:** Can't store 10,000 turning points in RAM.

**Solution:** Use delegates to stream turning points to external flash/EEPROM.

.. code-block:: c

   bool tp_set_flash(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s *tp)
   {
       uint32_t addr = FLASH_TP_BASE + (tp_pos - 1) * sizeof(*tp);
       flash_write(addr, tp, sizeof(*tp));
       return true;
   }

   bool tp_get_flash(rfc_ctx_s *ctx, size_t tp_pos, rfc_value_tuple_s **tp)
   {
       static rfc_value_tuple_s cache[4];  // Small cache
       uint32_t addr = FLASH_TP_BASE + (tp_pos - 1) * sizeof(*tp);
       flash_read(addr, &cache[tp_pos % 4], sizeof(*tp));
       *tp = &cache[tp_pos % 4];
       return true;
   }

2. Real-Time Logging
--------------------

**Problem:** Need to log every cycle as it's detected.

**Solution:** Use cycle find delegate to intercept counting.

.. code-block:: c

   void cycle_find_with_logging(rfc_ctx_s *ctx, rfc_flags_e flags)
   {
       // Call default implementation
       // (you'd need to store the original function pointer)
       default_cycle_find(ctx, flags);

       // Log after each cycle
       if (ctx->full_inc > prev_full_inc)
       {
           log_cycle(ctx->class_from, ctx->class_to);
       }
   }

3. Material Database Integration
---------------------------------

**Problem:** Damage depends on material properties from database.

**Solution:** Use damage calculation delegate with database lookup.

.. code-block:: c

   bool damage_from_db(rfc_ctx_s *ctx, unsigned from_class, unsigned to_class,
                      double *damage, double *Sa_ret)
   {
       // Custom function to get temperature from context, context holds time positional data for this purpose
       float temperature = get_current_temperature(ctx);
       material_db_s *db = get_database_for_temperature(temperature);

       double range = abs(to_class - from_class) * ctx->class_width;
       double Sa = range / 2.0;

       // Query database for Wöhler parameters at current temperature
       double k = db_query_woehler_k(db, temperature);
       double Sx = db_query_woehler_Sx(db, temperature);

       // Calculate damage
       double Nf = pow(Sa / Sx, -k);
       *damage = 1.0 / Nf;
       *Sa_ret = Sa;

       return true;
   }

Troubleshooting
===============

Common Issues
-------------

**1. Delegate Not Called**

   .. code-block:: text

      Problem: Set delegate but it's never called.

   **Solutions:**

   - Ensure ``RFC_USE_DELEGATES`` is enabled at compile time
   - For ``cycle_find_fcn``, set ``counting_method = RFC_COUNTING_METHOD_DELEGATED``
   - Check that feature flag is enabled (e.g., ``RFC_TP_SUPPORT`` for tp delegates)

**2. Segmentation Fault**

   .. code-block:: text

      Problem: Crash when delegate is called.

   **Solutions:**

   - Ensure delegate function signature matches typedef exactly
   - Check for NULL pointer dereferences in your delegate
   - Verify ``ctx->internal.obj`` contains valid data

**3. Infinite Recursion**

   .. code-block:: text

      Problem: Stack overflow in delegate.

   **Solutions:**

   - Don't call RFC functions that trigger the same delegate
   - If calling default implementation, store original function pointer separately

**4. Memory Leaks**

   .. code-block:: text

      Problem: Memory not freed from delegate allocations.

   **Solutions:**

   - Store allocations in a structure referenced by ``ctx->internal.obj``
   - Free in cleanup code after ``RFC_deinit()``

Limitations
===========

Current Limitations
-------------------

1. **No Memory Allocation Delegates**

   The library uses standard ``malloc()/free()``. Custom allocators require modifying source or using external wrapping.

2. **Single Context Per Delegate**

   Each ``rfc_ctx_s`` has its own delegates. Can't share one delegate across multiple contexts easily without storing context mapping.

3. **C Function Pointers Only**

   Delegates must be C-style function pointers. C++ functors/std::function not directly supported (use static wrappers).

4. **Not Thread-Safe**

   Delegates are called without synchronization. Multi-threaded usage requires external locking.

See Also
========

- `turning_points.rst <turning_points.rst>`_ - External turning point storage details
- `spread_damage.rst <spread_damage.rst>`_ - Custom damage spreading
- `cpp_wrapper.rst <cpp_wrapper.rst>`_ - C++ integration with delegates
- `minimal_build.rst <minimal_build.rst>`_ - Why delegates are disabled in RFC_MINIMAL
- `examples.rst <examples.rst>`_ - More practical examples

References
==========

The delegate pattern is a common technique in embedded systems and extensible libraries:

- **Design Patterns** - Strategy pattern (similar concept)
- **Embedded C** - Function pointer callbacks for hardware abstraction
- **Event-Driven Systems** - Observer pattern with function pointers

For delegate usage in the rainflow context, consult:

- ``rainflow.h`` - Delegate typedef definitions (lines 585-605)
- ``rainflow.c`` - Delegate invocation examples (search for ``RFC_USE_DELEGATES``)
- Test files - ``rfc_test.c`` for delegate testing examples
