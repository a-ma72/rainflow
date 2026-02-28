========
Features
========

The rainflow counting package provides a comprehensive set of features for
fatigue analysis. This page describes the key capabilities and optional
features available.

Modular Architecture
====================

Two-Layer Design
----------------

The package uses a modular two-layer architecture:

**Layer 1: Core C Module** (``rainflow.c`` / ``rainflow.h``)
   Contains all necessary functions for rainflow counting and histogram
   extraction. Features can be selectively enabled at compile time.

**Layer 2: C++ Wrapper** (``rainflow.hpp``)
   Encapsulates C functions in a namespace and provides a template class
   ``Rainflow`` for object-oriented access and inheritance. Includes
   container-based turning point storage.

Compile-Time Feature Selection
-------------------------------

You can customize the build by enabling specific features:

``RFC_MINIMAL``
   Core functions for rainflow counting only (suitable for µControllers)

``RFC_TP_SUPPORT``
   Turning point storage

``RFC_HCM_SUPPORT``
   HCM algorithm (Clormann/Seeger 3-point method)

``RFC_ASTM_SUPPORT``
   ASTM E 1049 (2011) algorithm variant (3-point method with specific residue handling)

``RFC_AT_SUPPORT``
   User-defined amplitude transformation (Haigh diagram)

``RFC_DH_SUPPORT``
   Damage history storage

``RFC_USE_DELEGATES``
   Delegates for core functions to implement custom behavior

``RFC_GLOBAL_EXTREMA``
   Track global data extrema during counting

``RFC_DAMAGE_FAST``
   Use lookup tables for damage and amplitude transformation

``RFC_EXPORT_MEX``
   Export mexFunction() for MATLAB integration

``RFC_EXPORT_PY``
   Export Python extension module

``RFC_UNIT_TEST``
   Build executable for unit testing

.. note::

   You can use COAN_ to tidy up code and remove unwanted features. The
   minimal version uses ``RFC_MINIMAL`` only.

.. _COAN: http://coan2.sourceforge.net/

Core Features
=============

Streaming Capability
--------------------

Process data flexibly:

- **At once** - Feed entire dataset in one call
- **In packages** - Process data in chunks
- **Sample-wise** - Feed one sample at a time

This allows integration with real-time systems and memory-constrained
environments.

Dynamic Class Management
------------------------

- Configurable class width (bin size)
- Dynamic class range expansion (``RFC_FLAGS_AUTORESIZE``)
- Automatic class count adjustment when data exceeds initial range
- Class offset for custom histogram boundaries

Counting Methods
----------------

Three standard counting methods are supported:

**4-Point Method (Default)**
   Standard rainflow algorithm as described in :doc:`algorithm`

**HCM (Clormann/Seeger)**
   3-point Hysteresis Counting Method based on material mechanics

**ASTM Method**
   Variant specified in ASTM E 1049 (2011)

Fatigue Analysis
================

Wöhler Curve (S-N Curve)
-------------------------

Configure material fatigue behavior with:

- Up to two slopes (``k`` and ``k2``)
- Fatigue limit (``sd``, ``nd``)
- Endurance limit (``sx``, ``nx``)
- Omission threshold (cycles below threshold ignored)

Damage Accumulation
-------------------

Four variants of Miner's rule:

**Elementary Miner**
   D = Σ(n_i / N_i)

**Original Miner**
   Considers load sequence effects

**Modified Miner**
   Adjusts for non-linear damage accumulation

**Consistent Miner**
   Ensures consistent damage calculation with material degradation

In-Time Analysis
----------------

Track damage accumulation over time:

- **Damage History** - Damage at each data point
- **Damage Indicator** - Real-time fatigue state (consistent Miner)
- **Turning Points** - Marked with assigned damage values

Result Extraction
=================

Histograms
----------

**Rainflow Matrix (RFM)**
   2D histogram of cycle amplitudes vs. mean values

**Level Crossing (LC)**
   1D histogram of stress/strain level crossings

**Range Pair (RP)**
   1D histogram of cycle ranges

Turning Points
--------------

Access complete turning point information:

- Position in input stream
- Value (stress/strain level)
- Associated damage
- Pair marking (for closed hysteresis)

Residue Processing
------------------

Multiple methods for handling unclosed cycles:

- **DIN 45667** - German standard method
- **ASTM halfcycle** - Count as 0.5 cycles
- **ASTM fullcycle** - Count as full cycles
- **Second run** - Re-feed residue
- **HCM** - Apply Clormann/Seeger method

Advanced Features
=================

Amplitude Transformation
------------------------

Support for mean stress effects (Haigh diagram):

- **FKM symmetrical** - Standard FKM guideline
- **FKM non-symmetrical** - Extended FKM method
- **User-defined** - Custom transformation functions

Lookup Tables
-------------

Accelerate damage calculation:

- Pre-computed damage values (``RFC_DAMAGE_FAST``)
- Amplitude transformation tables
- Automatic table generation

Damage History
--------------

Two storage modes:

**Compact History**
   Stores damage only at turning points with hysteresis pairs

**Uncompressed History**
   Full damage timeline matching input data length

Custom Delegates
----------------

Implement custom behavior via function pointers:

- Custom memory allocation
- User-defined counting logic
- Custom residue handling
- Extensible architecture

Data Conversions
================

Built-in conversion utilities:

- RFM → LC (Rainflow matrix to level crossing)
- RFM → RP (Rainflow matrix to range pairs)
- RFM → Damage (Direct damage from matrix)
- RP → Damage (Damage from range pairs, all Miner variants)

Language Bindings
=================

The package provides bindings for multiple environments:

**Python**
   Full-featured extension module (``rfcnt``)

   - NumPy array support
   - Pythonic API
   - Comprehensive error handling

**MATLAB**
   MEX interface for MATLAB integration

   - Native MATLAB array handling
   - Command-line and script usage

**C/C++**
   Direct library usage

   - Header-only C++ wrapper
   - C99-compliant core
   - Portable across platforms

See Also
========

- `installation.rst <installation.rst>`_ - Build configuration
- `examples.rst <examples.rst>`_ - Usage examples
- `algorithm.rst <algorithm.rst>`_ - Algorithm details
