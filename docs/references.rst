==========
References
==========

This page contains the bibliography and citations for the standards,
publications, and resources referenced in the rainflow counting package.

Standards
=========

.. [1] **ASTM E 1049 - Standard Practices for Cycle Counting in Fatigue Analysis**

   ASTM Standard E 1049, 1985 (2011).
   West Conshohocken, PA: ASTM International, 2011.

   This standard provides comprehensive guidance on cycle counting methods
   for fatigue analysis, including the 4-point rainflow algorithm
   implemented in this package.

Scientific Publications
=======================

.. [2] **Rainflow - HCM / Ein Hysteresisschleifen-Zählalgorithmus auf werkstoffmechanischer Grundlage**

   U.H. Clormann, T. Seeger

   1985, TU Darmstadt, Fachgebiet Werkstoffmechanik

   Introduces the Hysteresis Counting Method (HCM), a 3-point algorithm
   based on material mechanics principles. This method is implemented as
   an alternative counting approach in the package.

.. [3] **Zählverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen**

   FVA-Richtlinie, 2010.

   Available: https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf

   German Research Association for Drive Technology (FVA) guideline
   describing counting methods for creating collectives and matrices
   from time series data.

.. [4] **Rainflow Counting**

   Siemens Product Lifecycle Management Software Inc., 2018.

   Available: https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093

   Technical documentation on rainflow counting implementation and
   practical applications in fatigue testing.

.. [5] **Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation**

   G. Marsh

   International Journal of Fatigue 82 (2016) 757-765

   DOI: https://doi.org/10.1016/j.ijfatigue.2015.10.007

   Comprehensive review of different methods for handling the residue
   (unclosed cycles) in rainflow counting, comparing accuracy and
   computational efficiency.

Books and Monographs
====================

.. [6] **Betriebsfestigkeit - Verfahren und Daten zur Bauteilberechnung**

   Haibach, Erwin

   Springer Verlag

   Comprehensive German-language reference on operational strength
   (fatigue under service loading), covering calculation methods and
   material data for component design.

.. [7] **Schaedigungsbasierte Hysteresefilter**

   Hack, M.

   D386 (Diss. Univ. Kaiserslautern), Shaker Verlag Aachen, 1998

   ISBN 3-8265-3936-2

   Dissertation on damage-based hysteresis filters for fatigue analysis.

.. [8] **Hysteresis and Phase Transition**

   Brokate, M.; Sprekels, J.

   Applied Mathematical Sciences 121, Springer, New York, 1996

   Mathematical treatment of hysteresis phenomena with applications to
   material behavior and phase transitions.

.. [9] **Rainflow counting and energy dissipation in elastoplasticity**

   Brokate, M.; Dressler, K.; Krejci, P.

   Eur. J. Mech. A/Solids 15, pp. 705-737, 1996

   Theoretical foundation connecting rainflow counting to energy
   dissipation in elastoplastic materials.

.. [10] **Multivariate Density Estimation: Theory, Practice and Visualization**

   Scott, D.

   New York, Chichester, Wiley & Sons, 1992

   Statistical methods relevant to histogram-based analysis and
   density estimation in cycle counting.

.. [11] **Werkstoffmechanik - Bauteile sicher beurteilen und Werkstoffe richtig einsetzen**

   Ralf Bürgel, Hans Albert Richard, André Riemer

   Springer Fachmedien Wiesbaden 2005, 2014

   German textbook on materials mechanics, covering safe assessment
   of components and proper material selection.

.. [12] **Zählverfahren und Lastannahme in der Betriebsfestigkeit**

   Michael Köhler, Sven Jenne / Kurt Pötter, Harald Zenner

   Springer-Verlag Berlin Heidelberg 2012

   Comprehensive German reference on counting methods and load assumptions
   in operational fatigue strength analysis.

.. [13] **DIN 45667-1:2025-12** Mechanische Schwingungen - Zählverfahren für die Betriebsfestigkeit - Teil 1: Vorverarbeitung

   Beuth Verlag GmbH, Berlin

   Current standard defining signal preprocessing, specifically the extraction of
   turning points and discretization (binning) for fatigue analysis.

.. [14] **DIN 45667-2:2025-12** Mechanische Schwingungen - Zählverfahren für die Betriebsfestigkeit - Teil 2: Zählverfahren

   Beuth Verlag GmbH, Berlin

   Current standard defining cycle counting algorithms (Rainflow) and residue handling
   (closed vs. open cycles) for fatigue analysis.

Online Resources
================

COAN - C/C++ Configuration Preprocessor
----------------------------------------

http://coan2.sourceforge.net/

A source code preprocessing tool used to conditionally compile different
feature sets. Can be used to create minimal versions of the rainflow
library with only required features enabled.

GitHub Repository
-----------------

https://github.com/a-ma72/rainflow

Official repository containing the source code, issue tracker, and
release downloads.

Related Standards and Guidelines
================================

FKM Guideline
-------------

Forschungskuratorium Maschinenbau (FKM) - Analytical strength assessment
of components in mechanical engineering.

German standard for analytical strength assessment, including methods
for mean stress correction (Haigh diagram) implemented in the amplitude
transformation features.

Citation
========

If you use this software in your research, please cite:

.. code-block:: bibtex

   @software{rainflow,
     author = {Andreas Martin},
     title = {Rainflow Counting Algorithm},
     year = {2026},
     url = {https://github.com/a-ma72/rainflow},
     version = {0.5.2}
   }

Contributing
============

Contributions to improve documentation, add examples, or fix issues are
welcome. Please see the GitHub repository for contribution guidelines.

See Also
========

- `algorithm.rst <algorithm.rst>`_ - How these methods are implemented
- `features.rst <features.rst>`_ - Which standards/methods are supported
