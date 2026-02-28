==================================
Rainflow Documentation Index
==================================

Welcome to the comprehensive documentation for the Rainflow Counting library.

Getting Started
===============

New to Rainflow counting? Start here:

1. `algorithm.rst <algorithm.rst>`_ - Understanding the rainflow algorithm
2. `installation.rst <installation.rst>`_ - Building and installing
3. `examples.rst <examples.rst>`_ - Quick code examples

Core Documentation
==================

Algorithm & Methods
-------------------

`algorithm.rst <algorithm.rst>`_
   Core rainflow counting algorithm explanation with 4-point method

`astm_method.rst <astm_method.rst>`_
   ASTM E 1049 method vs. 4-point algorithm - performance comparison

`residue_methods.rst <residue_methods.rst>`_
   Methods for handling unclosed cycles (residue)

Implementation Guides
---------------------

`cpp_wrapper.rst <cpp_wrapper.rst>`_
   C++ wrapper guide - RAII, templates, STL integration

`minimal_build.rst <minimal_build.rst>`_
   Building minimal version for embedded systems (RFC_MINIMAL)

`delegates.rst <delegates.rst>`_
   Custom behavior via delegate functions

`turning_points.rst <turning_points.rst>`_
   External turning point storage and management

Advanced Features
-----------------

`damage_history.rst <damage_history.rst>`_
   Tracking damage accumulation over time

`spread_damage.rst <spread_damage.rst>`_
   Transient damage distribution methods

`tp_prune.rst <tp_prune.rst>`_
   RFC_tp_prune() logic and usage

Reference Documentation
=======================

`features.rst <features.rst>`_
   Complete feature list and compile-time options

`installation.rst <installation.rst>`_
   Platform-specific build instructions

`examples.rst <examples.rst>`_
   Code examples in Python, C/C++, and MATLAB

`references.rst <references.rst>`_
   Bibliography and standards citations

Quick Reference
===============

Compile-Time Options
--------------------

.. code-block:: bash

   # Full featured build
   cmake -S. -Bbuild \
       -DRFC_TP_SUPPORT=1 \
       -DRFC_HCM_SUPPORT=1 \
       -DRFC_ASTM_SUPPORT=1 \
       -DRFC_DH_SUPPORT=1 \
       -DRFC_USE_DELEGATES=1

   # Minimal build (embedded systems)
   cmake -S. -Bbuild -DRFC_MINIMAL=1

Common Tasks
------------

**Enable ASTM counting:**
   See `astm_method.rst <astm_method.rst>`_

**Track damage over time:**
   See `damage_history.rst <damage_history.rst>`_

**Handle residue:**
   See `residue_methods.rst <residue_methods.rst>`_

**Custom memory allocation:**
   See `delegates.rst <delegates.rst>`_

**Prune turning points:**
   See `tp_prune.rst <tp_prune.rst>`_

Documentation by Topic
======================

By Use Case
-----------

**Embedded Systems**
   - `minimal_build.rst <minimal_build.rst>`_ - Minimal configuration
   - `delegates.rst <delegates.rst>`_ - Custom memory management
   - `turning_points.rst <turning_points.rst>`_ - External storage

**Performance Optimization**
   - `astm_method.rst <astm_method.rst>`_ - Fastest counting method
   - `tp_prune.rst <tp_prune.rst>`_ - Memory management
   - `features.rst <features.rst>`_ - Feature selection

**Research & Analysis**
   - `damage_history.rst <damage_history.rst>`_ - Time-based damage
   - `spread_damage.rst <spread_damage.rst>`_ - Damage distribution
   - `residue_methods.rst <residue_methods.rst>`_ - Residue handling

**Integration**
   - `cpp_wrapper.rst <cpp_wrapper.rst>`_ - C++ projects
   - `examples.rst <examples.rst>`_ - Python, MATLAB
   - `delegates.rst <delegates.rst>`_ - Custom callbacks

By Programming Language
-----------------------

**C**
   - `algorithm.rst <algorithm.rst>`_ - Core concepts
   - `delegates.rst <delegates.rst>`_ - Function pointers
   - `minimal_build.rst <minimal_build.rst>`_ - Minimal API

**C++**
   - `cpp_wrapper.rst <cpp_wrapper.rst>`_ - Complete C++ guide
   - `turning_points.rst <turning_points.rst>`_ - STL containers
   - `examples.rst <examples.rst>`_ - C++ examples

**Python**
   - `examples.rst <examples.rst>`_ - Python API
   - `damage_history.rst <damage_history.rst>`_ - Visualization
   - `installation.rst <installation.rst>`_ - pip install

By Feature
----------

**RFC_MINIMAL**
   `minimal_build.rst <minimal_build.rst>`_

**RFC_TP_SUPPORT**
   `turning_points.rst <turning_points.rst>`_

**RFC_HCM_SUPPORT**
   `algorithm.rst <algorithm.rst>`_ (HCM section)

**RFC_ASTM_SUPPORT**
   `astm_method.rst <astm_method.rst>`_

**RFC_DH_SUPPORT**
   `damage_history.rst <damage_history.rst>`_

**RFC_USE_DELEGATES**
   `delegates.rst <delegates.rst>`_

Glossary
========

**Rainflow Counting**
   Algorithm for extracting fatigue cycles from variable amplitude loading

**Turning Point**
   Local maximum or minimum in the load signal

**Residue**
   Unclosed cycles remaining after rainflow counting

**Damage History**
   Cumulative damage tracked at each time point

**Wöhler Curve (S-N Curve)**
   Material fatigue properties relating stress to cycles-to-failure

**Miner's Rule**
   Linear damage accumulation hypothesis

**HCM**
   Hysteresis Counting Method (Clormann/Seeger)

**ASTM E 1049**
   Standard for cycle counting in fatigue analysis

**Delegate**
   User-provided callback function for custom behavior

Contributing
============

Found an error or want to improve documentation?

- Report issues: https://github.com/a-ma72/rainflow/issues
- Submit PRs: https://github.com/a-ma72/rainflow/pulls

License
=======

See the main repository for license information.

---

**Last Updated:** 2026-01
**Documentation Version:** 0.5.2
