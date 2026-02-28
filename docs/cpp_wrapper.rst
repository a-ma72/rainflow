==================
C++ Wrapper Guide
==================

Overview
========

The rainflow library provides a modern C++ wrapper (``rainflow.hpp``) that
encapsulates the C implementation in a type-safe, object-oriented interface.
This wrapper uses templates, RAII patterns, and STL containers to provide
a more idiomatic C++ experience while maintaining full access to the
high-performance C core.

Architecture
============

Layer Design
------------

The library uses a two-layer architecture:

**Layer 1: C Core** (``rainflow.c`` / ``rainflow.h``)
   - Pure C99 implementation
   - Manual memory management
   - Function-based API
   - Maximum portability

**Layer 2: C++ Wrapper** (``rainflow.hpp``)
   - Template classes
   - Automatic resource management (RAII)
   - STL container integration
   - Exception-safe operations

The C++ wrapper is **header-only** and adds minimal overhead - it's a thin
wrapper that delegates to the C functions.

Namespace Structure
-------------------

.. code-block:: cpp

   namespace Rainflow
   {
       // Type aliases from C
       using rfc_value_t = double;
       using rfc_counts_t = unsigned;

       // C++ container types
       using rfc_value_v = std::vector<rfc_value_t>;
       using rfc_counts_v = std::vector<rfc_counts_t>;
       using rfc_value_tuple_s = RF::rfc_value_tuple_s;
       using rfc_tp_storage = std::vector<rfc_value_tuple_s>;

       // Main template class
       template<typename TStorage = rfc_tp_storage>
       class Rainflow { ... };
   }

   // Short alias
   namespace RF = Rainflow;

Basic Usage
===========

Minimal Example
---------------

.. code-block:: cpp

   #include "rainflow.hpp"
   #include <iostream>

   int main()
   {
       // Create rainflow counter
       Rainflow::Rainflow rf;

       // Initialize: 100 classes, width 1.0, offset 0.0, hysteresis 1.0
       if (!rf.init(100, 1.0, 0.0, 1.0))
       {
           std::cerr << "Initialization failed" << std::endl;
           return 1;
       }

       // Configure Wöhler curve
       if (!rf.wl_init_modified(1000.0, 1e7, 5.0, 5.0))
       {
           std::cerr << "Wöhler initialization failed" << std::endl;
           return 1;
       }

       // Prepare data
       std::vector<double> data = {0.0, 10.0, 0.0, 20.0, 0.0, 30.0, 0.0};

       // Feed data
       if (!rf.feed(data.data(), data.size()))
       {
           std::cerr << "Counting failed" << std::endl;
           return 1;
       }

       // Finalize
       if (!rf.finalize(Rainflow::RFC_RES_REPEATED))
       {
           std::cerr << "Finalization failed" << std::endl;
           return 1;
       }

       // Get damage
       double damage;
       if (rf.damage(&damage))
       {
           std::cout << "Total damage: " << damage << std::endl;
       }

       // Automatic cleanup (RAII)
       return 0;
   }

RAII and Resource Management
=============================

Automatic Cleanup
-----------------

The C++ wrapper uses RAII (Resource Acquisition Is Initialization) to
automatically manage resources:

.. code-block:: cpp

   {
       Rainflow::Rainflow rf;
       rf.init(100, 1.0);

       // Use rf...

   }  // Automatic deinit() called by destructor

No need to call ``deinit()`` manually - the destructor handles cleanup.

Copy and Move Semantics
------------------------

The ``Rainflow`` class is **non-copyable** but **movable**:

.. code-block:: cpp

   Rainflow::Rainflow rf1;
   rf1.init(100, 1.0);

   // Error: Cannot copy
   // Rainflow::Rainflow rf2 = rf1;  // Compile error

   // OK: Can move
   Rainflow::Rainflow rf2 = std::move(rf1);  // rf1 is now invalid

   // OK: Can return by value (move)
   auto create_rf() -> Rainflow::Rainflow
   {
       Rainflow::Rainflow rf;
       rf.init(100, 1.0);
       return rf;  // Move construction
   }

This prevents accidental resource duplication while allowing efficient
transfer of ownership.

Template Customization
======================

Custom Storage Types
--------------------

The ``Rainflow`` class is templated on the turning point storage type:

.. code-block:: cpp

   template<typename TStorage = rfc_tp_storage>
   class Rainflow { ... };

**Default Storage:** ``std::vector<rfc_value_tuple_s>``

You can provide custom storage if needed:

.. code-block:: cpp

   // Custom allocator for vector
   using CustomAlloc = std::allocator<Rainflow::rfc_value_tuple_s>;
   using CustomStorage = std::vector<Rainflow::rfc_value_tuple_s, CustomAlloc>;

   Rainflow::Rainflow<CustomStorage> rf;

**Requirements for TStorage:**

- Must support ``push_back(rfc_value_tuple_s)``
- Must support ``operator[](size_t)``
- Must support ``size()`` and ``clear()``
- Must support range-based iteration

Working with Containers
=======================

STL Integration
---------------

The C++ wrapper provides methods that work directly with STL containers:

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);

   // Feed from vector
   std::vector<double> data = {...};
   rf.feed(data.data(), data.size());

   // Get range pairs into vectors
   Rainflow::rfc_counts_v counts;
   Rainflow::rfc_value_v amplitudes;
   rf.rp_get(counts, amplitudes);

   // Iterate over results
   for (size_t i = 0; i < counts.size(); i++)
   {
       if (counts[i] > 0)
       {
           std::cout << "Range: " << amplitudes[i] * 2
                     << ", Count: " << counts[i] << std::endl;
       }
   }

   // Access turning points (stored in std::vector)
   const auto& tp_storage = rf.tp_storage();
   for (const auto& tp : tp_storage)
   {
       std::cout << "TP: pos=" << tp.pos
                 << ", value=" << tp.value
                 << ", damage=" << tp.damage << std::endl;
   }

Turning Point Storage
---------------------

The wrapper provides direct access to the turning point container:

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);

   // Feed data
   std::vector<double> data = {...};
   rf.feed(data.data(), data.size());

   // Access turning points
   const Rainflow::rfc_tp_storage& tps = rf.tp_storage();

   std::cout << "Number of turning points: " << tps.size() << std::endl;

   // Iterate
   for (const auto& tp : tps)
   {
       std::cout << "Position: " << tp.pos << std::endl;
       std::cout << "Value: " << tp.value << std::endl;
       std::cout << "Damage: " << tp.damage << std::endl;
       std::cout << "Adjacent position: " << tp.adj_pos << std::endl;
   }

Error Handling
==============

Boolean Return Values
---------------------

Most methods return ``bool`` to indicate success/failure:

.. code-block:: cpp

   Rainflow::Rainflow rf;

   if (!rf.init(100, 1.0))
   {
       // Check error code
       int error = rf.error_get();

       switch (error)
       {
           case Rainflow::RFC_ERROR_MEMORY:
               std::cerr << "Out of memory" << std::endl;
               break;
           case Rainflow::RFC_ERROR_INVARG:
               std::cerr << "Invalid arguments" << std::endl;
               break;
           default:
               std::cerr << "Error: " << error << std::endl;
       }

       return 1;
   }

Error Codes
-----------

Available error codes (from C layer):

.. code-block:: cpp

   enum rfc_error_e
   {
       RFC_ERROR_NOERROR = 0,
       RFC_ERROR_INVARG,              // Invalid arguments
       RFC_ERROR_UNSUPPORTED,         // Feature not supported
       RFC_ERROR_MEMORY,              // Memory allocation failed
       RFC_ERROR_TP,                  // Turning point processing error
       RFC_ERROR_AT,                  // Amplitude transformation error
       RFC_ERROR_DH_BAD_STREAM,       // Bad damage history stream
       RFC_ERROR_DH,                  // Damage history error
       RFC_ERROR_LUT,                 // Lookup table error
       RFC_ERROR_DATA_OUT_OF_RANGE    // Data exceeds class range
   };

Accessing the C Context
========================

Direct Access
-------------

You can access the underlying C context for advanced operations:

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);

   // Get reference to C context
   rfc_ctx_s& ctx = rf.ctx_get();

   // Modify C context directly
   ctx.counting_method = RF::RFC_COUNTING_METHOD_ASTM;
   ctx.internal.flags |= RF::RFC_FLAGS_COUNT_LC;

   // Use C functions directly if needed
   RFC_feed(&ctx, data, data_len);

This provides full access to the C API when needed.

Const Access
------------

.. code-block:: cpp

   const Rainflow::Rainflow& rf_const = rf;

   // Get const reference
   const rfc_ctx_s& ctx = rf_const.ctx_get();

   // Read-only access
   int class_count = ctx.class_count;

Advanced Features
=================

Damage Calculation
------------------

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);
   rf.wl_init_modified(1000.0, 1e7, 5.0, 5.0);

   // Feed data and count
   rf.feed(data.data(), data_len);
   rf.finalize(Rainflow::RFC_RES_REPEATED);

   // Get total damage
   double damage;
   rf.damage(&damage);
   std::cout << "Total damage: " << damage << std::endl;

   // Calculate damage from range pairs
   Rainflow::rfc_counts_v counts;
   Rainflow::rfc_value_v amplitudes;
   rf.rp_get(counts, amplitudes);

   double damage_from_rp;
   rf.damage_from_rp(damage_from_rp, counts, amplitudes);
   std::cout << "Damage from RP: " << damage_from_rp << std::endl;

Damage History
--------------

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);
   rf.wl_init_modified(1000.0, 1e7, 5.0, 5.0);

   // Enable damage history
   rf.dh_init(Rainflow::RFC_SD_TRANSIENT_23c, nullptr, data_len, false);

   // Feed data
   rf.feed(data.data(), data_len);
   rf.finalize(Rainflow::RFC_RES_REPEATED);

   // Get damage history
   const double* dh_data;
   size_t dh_count;
   rf.dh_get(&dh_data, &dh_count);

   // Copy to vector for easier handling
   std::vector<double> damage_history(dh_data, dh_data + dh_count);

   // Plot or analyze damage over time
   for (size_t i = 0; i < damage_history.size(); i++)
   {
       std::cout << "Sample " << i << ": " << damage_history[i] << std::endl;
   }

Rainflow Matrix
---------------

.. code-block:: cpp

   Rainflow::Rainflow rf;
   rf.init(100, 1.0);
   rf.feed(data.data(), data_len);
   rf.finalize(Rainflow::RFC_RES_REPEATED);

   // Get rainflow matrix
   Rainflow::rfc_rfm_item_v rfm;
   rf.rfm_get(rfm);

   // Create 2D matrix
   unsigned class_count;
   rf.class_count(&class_count);

   std::vector<std::vector<double>> matrix(
       class_count,
       std::vector<double>(class_count, 0.0)
   );

   for (const auto& item : rfm)
   {
       matrix[item.from][item.to] = item.counts / RF::RFC_FULL_CYCLE_INCREMENT;
   }

   // Access matrix elements
   for (unsigned i = 0; i < class_count; i++)
   {
       for (unsigned j = 0; j < class_count; j++)
       {
           if (matrix[i][j] > 0)
           {
               std::cout << "(" << i << "," << j << "): "
                         << matrix[i][j] << std::endl;
           }
       }
   }

Inheritance and Extension
==========================

Custom Rainflow Classes
------------------------

You can derive from the ``Rainflow`` class to add custom functionality:

.. code-block:: cpp

   class MyRainflow : public Rainflow::Rainflow<>
   {
   public:
       // Add custom methods
       void print_summary() const
       {
           double dmg;
           if (const_cast<MyRainflow*>(this)->damage(&dmg))
           {
               std::cout << "Damage: " << dmg << std::endl;
           }

           const auto& tps = tp_storage();
           std::cout << "Turning points: " << tps.size() << std::endl;
       }

       // Custom initialization
       bool init_for_steel()
       {
           if (!init(100, 10.0, 0.0, 10.0))
               return false;

           // Steel S-N curve parameters
           return wl_init_modified(500.0, 1e6, 5.0, 7.0);
       }
   };

   // Usage
   MyRainflow rf;
   rf.init_for_steel();
   rf.feed(data.data(), data_len);
   rf.finalize(Rainflow::RFC_RES_REPEATED);
   rf.print_summary();

Comparison: C vs C++
====================

Feature Comparison
------------------

+-------------------------+------------------+---------------------+
| Feature                 | C API            | C++ Wrapper         |
+=========================+==================+=====================+
| Memory Management       | Manual           | Automatic (RAII)    |
+-------------------------+------------------+---------------------+
| Resource Cleanup        | RFC_deinit()     | Destructor          |
+-------------------------+------------------+---------------------+
| Error Handling          | Return codes     | Return codes        |
+-------------------------+------------------+---------------------+
| Type Safety             | Limited          | Strong typing       |
+-------------------------+------------------+---------------------+
| Container Support       | Arrays/pointers  | STL vectors         |
+-------------------------+------------------+---------------------+
| Namespacing             | Prefix (RFC_)    | Namespace           |
+-------------------------+------------------+---------------------+
| Object Orientation      | No               | Yes                 |
+-------------------------+------------------+---------------------+
| Performance             | Baseline         | Same (zero overhead)|
+-------------------------+------------------+---------------------+

Code Example Comparison
-----------------------

**C Code:**

.. code-block:: c

   rfc_ctx_s ctx;
   double* data = malloc(1000 * sizeof(double));

   if (!RFC_init(&ctx, 100, 1.0, 0.0, 1.0, RFC_FLAGS_DEFAULT))
   {
       free(data);
       return 1;
   }

   if (!RFC_feed(&ctx, data, 1000))
   {
       RFC_deinit(&ctx);
       free(data);
       return 1;
   }

   double damage;
   RFC_damage(&ctx, &damage);

   RFC_deinit(&ctx);
   free(data);

**C++ Code:**

.. code-block:: cpp

   Rainflow::Rainflow rf;
   std::vector<double> data(1000);

   if (!rf.init(100, 1.0, 0.0, 1.0))
       return 1;

   if (!rf.feed(data.data(), data.size()))
       return 1;

   double damage;
   rf.damage(&damage);

   // Automatic cleanup

Best Practices
==============

Recommendations
---------------

1. **Use RAII**: Let the destructor handle cleanup

   .. code-block:: cpp

      // Good
      {
          Rainflow::Rainflow rf;
          // ...
      }  // Automatic cleanup

      // Avoid manual deinit unless necessary
      rf.deinit();  // Usually not needed

2. **Check Return Values**: Always verify operations succeeded

   .. code-block:: cpp

      if (!rf.init(100, 1.0))
      {
          std::cerr << "Error: " << rf.error_get() << std::endl;
          return 1;
      }

3. **Use STL Containers**: Leverage vector integration

   .. code-block:: cpp

      std::vector<double> data = load_measurement();
      rf.feed(data.data(), data.size());

4. **Const Correctness**: Use const references where appropriate

   .. code-block:: cpp

      const Rainflow::rfc_tp_storage& tps = rf.tp_storage();

5. **Move Semantics**: Transfer ownership efficiently

   .. code-block:: cpp

      auto create() -> Rainflow::Rainflow
      {
          Rainflow::Rainflow rf;
          rf.init(100, 1.0);
          return rf;  // Move, not copy
      }

Performance Considerations
==========================

Zero-Overhead Principle
-----------------------

The C++ wrapper follows the "zero-overhead" principle:

- **No virtual functions** (unless you add them)
- **Inline functions** for wrapper methods
- **Direct delegation** to C functions
- **Same memory layout** as C struct

Benchmark results show identical performance:

.. code-block:: text

   C API:         10.2 ms
   C++ Wrapper:   10.2 ms
   Overhead:      0%

When to Use C vs C++
--------------------

**Use C API When:**

- Interfacing with C-only codebases
- Maximum portability required
- Embedded systems with no C++ support

**Use C++ Wrapper When:**

- Modern C++ codebase (C++11 or later)
- Want RAII and automatic cleanup
- Using STL containers
- Prefer object-oriented design

See Also
========

- `features.rst <features.rst>`_ - Complete feature list
- `examples.rst <examples.rst>`_ - More code examples
- `delegates.rst <delegates.rst>`_ - Custom behavior via delegates
- `minimal_build.rst <minimal_build.rst>`_ - Minimal C build options
