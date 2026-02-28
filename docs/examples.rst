========
Examples
========

This page provides practical examples of using the rainflow counting package
in different programming environments.

Python Examples
===============

Basic Rainflow Counting
------------------------

Simple example with synthetic load data:

.. code-block:: python

   import numpy as np
   import rfcnt

   # Generate synthetic load data
   data = np.array([0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0])

   # Perform rainflow counting
   result = rfcnt.rfc(
       data,
       class_width=0.5,
       class_count=10
   )

   # Access results
   print(f"Total damage: {result['damage']}")
   print(f"Range pairs:\n{result['rp']}")
   print(f"Level crossings:\n{result['lc']}")

With Wöhler Curve
-----------------

Configure fatigue parameters using a Wöhler (S-N) curve:

.. code-block:: python

   import numpy as np
   import rfcnt

   # Load measurement data
   data = np.loadtxt('strain_measurement.txt')

   # Define Wöhler curve parameters
   wl = {
       'sd': 0.0,      # Fatigue strength (0 = no limit)
       'nd': np.inf,   # Cycles at sd
       'sx': 1000.0,   # Point on sn curve (dividing the curve into k and k2)
       'nx': 1e7,      # Cycles at sx
       'k': 5.0,       # Slope above sx
       'k2': 7.0,      # Slope below sx (for two-slope curve)
       'omission': 50.0  # Ignore cycles below this value
   }

   # Perform counting with damage calculation
   result = rfcnt.rfc(
       data,
       class_width=10.0,
       wl=wl,
       auto_resize=True
   )

   print(f"Accumulated damage: {result['damage']:.6f}")
   print(f"Turning points: {len(result['tp'])}")

Advanced Options
----------------

Using different counting methods and options:

.. code-block:: python

   import rfcnt

   # HCM (Clormann/Seeger) method with damage history
   result_hcm = rfcnt.rfc(
       data,
       class_width=5.0,
       use_HCM=True,
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
       hysteresis=2.5
   )

   # ASTM method with specific residue handling
   result_astm = rfcnt.rfc(
       data,
       class_width=5.0,
       use_ASTM=True,
       residual_method=rfcnt.ResidualMethod.ASTM_FULLCYCLES
   )

   # With damage history
   result_dh = rfcnt.rfc(
       data,
       class_width=5.0,
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
       wl=wl
   )

   # Plot damage over time
   import matplotlib.pyplot as plt
   plt.plot(result_dh['dh'])
   plt.xlabel('Sample index')
   plt.ylabel('Accumulated damage')
   plt.title('Damage History')
   plt.show()

Damage from Range Pairs
------------------------

Calculate damage directly from existing range pair data:

.. code-block:: python

   import numpy as np
   import rfcnt

   # Existing range pair histogram
   # Sa: stress amplitudes, counts: cycle counts
   Sa = np.array([100, 150, 200, 250, 300])
   counts = np.array([1000, 500, 200, 50, 10])

   # Calculate damage
   damage = rfcnt.damage_from_rp(
       Sa,
       counts,
       wl={
           'sx': 1000,
           'nx': 1e7,
           'k': 5
       },
       method=rfcnt.RPDamageCalcMethod.MINER_ELEMENTARY
   )

   print(f"Total damage: {damage:.6e}")

Streaming/Chunked Processing
-----------------------------

Process large datasets in chunks to manage memory:

.. code-block:: python

   import rfcnt
   import numpy as np

   # Initialize counting context
   class_width = 10.0
   class_count = 100

   # Process first chunk
   chunk1 = np.loadtxt('data_part1.txt')
   result = rfcnt.rfc(chunk1, class_width=class_width, class_count=class_count)

   # Continue with more data
   # Note: For true streaming, you would need to preserve context
   # This is a simplified example

C/C++ Examples
==============

Basic C Usage
-------------

.. code-block:: c

   #include "rainflow.h"
   #include <stdio.h>

   int main(void)
   {
       rfc_ctx_s ctx;
       double data[] = {0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0};
       size_t data_len = sizeof(data) / sizeof(data[0]);
       double damage;

       // Initialize context
       if (!RFC_init(&ctx, 100, 0.5, 0.0, 0.5, RFC_FLAGS_DEFAULT))
       {
           fprintf(stderr, "Initialization failed\n");
           return 1;
       }

       // Configure Wöhler curve
       if (!RFC_wl_init_modified(&ctx, 1000.0, 1e7, 5.0, 5.0))
       {
           fprintf(stderr, "Wöhler init failed\n");
           RFC_deinit(&ctx);
           return 1;
       }

       // Feed data
       if (!RFC_feed(&ctx, data, data_len))
       {
           fprintf(stderr, "Counting failed\n");
           RFC_deinit(&ctx);
           return 1;
       }

       // Finalize and get results
       if (!RFC_finalize(&ctx, RFC_RES_REPEATED))
       {
           fprintf(stderr, "Finalization failed\n");
           RFC_deinit(&ctx);
           return 1;
       }

       // Get damage
       if (RFC_damage(&ctx, &damage))
       {
           printf("Total damage: %.6e\n", damage);
       }

       // Cleanup
       RFC_deinit(&ctx);
       return 0;
   }

C++ Wrapper Usage
-----------------

.. code-block:: cpp

   #include "rainflow.hpp"
   #include <iostream>
   #include <vector>

   int main()
   {
       // Create Rainflow object
       Rainflow rf;

       // Initialize with parameters
       if (!rf.init(100, 0.5, 0.0, 0.5))
       {
           std::cerr << "Initialization failed: "
                     << rf.error_get() << std::endl;
           return 1;
       }

       // Configure Wöhler curve
       if (!rf.wl_init_modified(1000.0, 1e7, 5.0, 5.0))
       {
           std::cerr << "Wöhler init failed" << std::endl;
           return 1;
       }

       // Prepare data
       std::vector<double> data = {0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0};

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

       // Get results
       double damage;
       if (rf.damage(&damage))
       {
           std::cout << "Total damage: " << damage << std::endl;
       }

       // Get range pairs
       Rainflow::rfc_counts_v counts;
       Rainflow::rfc_value_v amplitudes;
       if (rf.rp_get(counts, amplitudes))
       {
           std::cout << "Range pairs:" << std::endl;
           for (size_t i = 0; i < counts.size(); i++)
           {
               if (counts[i] > 0)
               {
                   std::cout << "  Range: " << amplitudes[i] * 2
                            << ", Count: " << counts[i] << std::endl;
               }
           }
       }

       rf.deinit();
       return 0;
   }

MATLAB Examples
===============

Basic MATLAB Usage
------------------

.. code-block:: matlab

   % Generate test data
   data = [0, 1, 0, 2, 0, 3, 0];

   % Perform rainflow counting
   result = rfc_mex(data, 'class_width', 0.5, 'class_count', 10);

   % Display results
   fprintf('Total damage: %.6e\n', result.damage);
   disp('Range pairs:');
   disp(result.rp);

With Wöhler Curve
-----------------

.. code-block:: matlab

   % Load measurement data
   data = load('strain_data.mat');

   % Define Wöhler parameters
   wl = struct('sx', 1000, 'nx', 1e7, 'k', 5, 'k2', 7, 'omission', 50);

   % Perform counting
   result = rfc_mex(data.strain, ...
                    'class_width', 10, ...
                    'wl', wl, ...
                    'auto_resize', true);

   fprintf('Damage: %.6f\n', result.damage);

   % Plot rainflow matrix
   figure;
   imagesc(result.rfm);
   colorbar;
   title('Rainflow Matrix');
   xlabel('To class');
   ylabel('From class');

Common Workflows
================

Complete Fatigue Analysis
--------------------------

End-to-end example of fatigue analysis workflow:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import rfcnt

   # 1. Load or generate data
   time = np.linspace(0, 100, 10000)
   strain = 100 * np.sin(2 * np.pi * 0.1 * time) + \
            50 * np.sin(2 * np.pi * 0.5 * time) + \
            np.random.normal(0, 10, len(time))

   # 2. Define material parameters
   wl = {
       'sx': 500,    # Stress amplitude at knee point
       'nx': 1e6,    # Cycles at knee point
       'k': 5,       # Slope above knee
       'k2': 9,      # Slope below knee
       'omission': 20  # Ignore small cycles
   }

   # 3. Perform rainflow counting
   result = rfcnt.rfc(
       strain,
       class_width=5.0,
       class_count=200,
       hysteresis=5.0,
       wl=wl,
       spread_damage=rfcnt.SDMethod.TRANSIENT_23c,
       auto_resize=True
   )

   # 4. Analyze results
   print(f"Total cycles: {np.sum(result['rp'][:, 1])}")
   print(f"Total damage: {result['damage']:.6e}")
   print(f"Remaining life: {1.0 / result['damage']:.1f} repetitions")

   # 5. Visualize
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))

   # Original signal with turning points
   axes[0, 0].plot(time, strain, 'b-', alpha=0.5, label='Signal')
   tp = result['tp']
   axes[0, 0].plot(tp[:, 0], tp[:, 1], 'ro', markersize=3, label='Turning points')
   axes[0, 0].set_xlabel('Time')
   axes[0, 0].set_ylabel('Strain')
   axes[0, 0].set_title('Input Signal')
   axes[0, 0].legend()
   axes[0, 0].grid(True)

   # Range pair histogram
   rp = result['rp']
   axes[0, 1].bar(rp[:, 0], rp[:, 1], width=5)
   axes[0, 1].set_xlabel('Range')
   axes[0, 1].set_ylabel('Count')
   axes[0, 1].set_title('Range Pair Histogram')
   axes[0, 1].grid(True)

   # Damage history
   axes[1, 0].plot(result['dh'])
   axes[1, 0].set_xlabel('Sample')
   axes[1, 0].set_ylabel('Accumulated Damage')
   axes[1, 0].set_title('Damage History')
   axes[1, 0].grid(True)

   # Rainflow matrix
   im = axes[1, 1].imshow(result['rfm'], origin='lower', aspect='auto')
   axes[1, 1].set_xlabel('To Class')
   axes[1, 1].set_ylabel('From Class')
   axes[1, 1].set_title('Rainflow Matrix')
   plt.colorbar(im, ax=axes[1, 1])

   plt.tight_layout()
   plt.savefig('fatigue_analysis.png', dpi=150)
   plt.show()

Running Built-in Examples
==========================

The Python package includes several example scripts:

.. code-block:: bash

   # Run all examples
   python -m rfcnt.run_examples

   # Run specific example
   python -m rfcnt.tests.examples.example_1

See Also
========

- `features.rst <features.rst>`_ - Detailed feature descriptions
- `algorithm.rst <algorithm.rst>`_ - Algorithm explanations
- `installation.rst <installation.rst>`_ - Setup instructions
