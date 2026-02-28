=====================
Rainflow Algorithm
=====================

Overview
========

The Rainflow counting algorithm is a widely used method in fatigue analysis
for extracting cycle information from time-dependent stress-strain curves.
This implementation uses the 4-point method in accordance with ASTM E 1049 and
DIN 45667.

Four Main Steps
===============

Rainflow Counting consists of four main steps:

1. **Hysteresis Filtering**

   Removes small oscillations below a specified threshold to focus on
   significant load variations.

2. **Peak-Valley Filtering**

   Identifies local maxima (peaks) and minima (valleys) in the filtered
   signal, representing turning points in the load history.

3. **Discretization**

   Maps continuous stress/strain values into discrete classes (bins) for
   histogram-based counting.

4. **Four Point Counting Method**

   The core algorithm that identifies closed cycles using a four-point pattern.

Four-Point Method
=================

The 4-point counting method identifies closed hysteresis loops by examining
sequences of four consecutive turning points::

                     * D
                    / \
             B *<--/
              / \ /
             /   * C
          \ /
           * A

A cycle B-C is considered **closed** if::

    min(B,C) >= min(A,D) && max(B,C) <= max(A,D)

When a cycle is closed:

- The slope B-C is counted as a full cycle
- Points B and C are removed from the residue
- The algorithm continues with the remaining points

This process continues until no more closed cycles can be found. The
remaining turning points form the **residue** (unclosed cycles).

Standards and References
========================

This implementation is fully documented in the following standards:

- **ASTM E 1049** - "Standard Practices for Cycle Counting in Fatigue Analysis" [1]_
- **HCM Method** (3-point) - Proposed by Clormann/Seeger [2]_
- **FVA Guidelines** - German Research Association for Drive Technology [3]_
- **DIN 45667-1:2025-12** - "Fatigue analysis of structures - Part 1: Rainflow counting method" [13]_
- **DIN 45667-2:2025-12** - "Fatigue analysis of structures - Part 2: Cycle counting method" [14]_

The implementation supports multiple counting methods:

- **4-point method** (default) - As described above
- **HCM (Hysteresis Counting Method)** - 3-point algorithm by Clormann/Seeger
- **ASTM method** - Variant specified in ASTM E 1049 (2011)

Residue Handling
================

After the main counting process, unclosed cycles (residue) can be processed
using various methods:

- **Ignore** - Discard residue (conservative approach)
- **DIN 45667** - According to German standard
- **ASTM halfcycle/fullcycle** - Half or full cycle counting, according to ASTM E 1049
- **Repeated** - Re-feed residue and count closed cycles
- **HCM** - Apply Clormann/Seeger method to residue

See `features.rst <features.rst>`_ for more details on residue handling options.

Differences between ASTM E 1049 vs. DIN 45667-2
------------------------------------------------

While the core Rainflow algorithm identifies closed hysteresis loops identically in both
standards, a critical difference exists in the handling of the residue
(the unclosed signal stream remaining at the end of the counting process).

**ASTM E 1049 (2011)** adopts a pragmatic approach, typically classifying these open loops as
"half cycles" (count = 0.5). While computationally straightforward, this method can
lead to underestimate fatigue damage, particularly in short signals representing
repetitive events (e.g., a single test track lap), as the largest stress range often
resides within the residue.

**DIN 45667-2:2025-12** adheres to a stricter physical definition. It mandates that all
significant stress ranges form complete hysteresis loops.
Therefore, it requires resolving the residue into full integer cycles (count = 1.0).
This is algorithmically achieved by signal reordering (starting the analysis at the absolute
global extremum) or block repetition (virtually appending the signal to itself).

Recommendation:

Select ASTM for strict compliance with US standards where half-cycles are explicitly permitted.

Select DIN (or enable "Close Residue") for conservative fatigue life prediction of repeating
load blocks, ensuring the global maximum stress range is fully weighted as a complete damage event.

Algorithm Complexity
====================

The algorithm operates in streaming mode, processing data points as they
arrive:

- **Time Complexity**: O(n) where n is the number of data points
- **Space Complexity**: O(r) where r is the residue size
- **Memory**: Configurable, supports both static and dynamic allocation

For large datasets, the implementation supports:

- Sample-by-sample processing
- Chunk-based streaming
- Dynamic class range expansion (auto-resize)

See Also
========

- `features.rst <features.rst>`_ - Implementation-specific features
- `references.rst <references.rst>`_ - Detailed citations
- `examples.rst <examples.rst>`_ - Practical usage examples
