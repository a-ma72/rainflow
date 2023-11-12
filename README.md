## Rainflow Counting Algorithm (4-point-methods), C99 compliant

### "Rainflow Counting" consists of four main steps:

  1. Hysteresis Filtering
  2. Peak-Valley Filtering
  3. Discretization
  4. Four Point Counting Method:

                     * D
                    / \       Closed, if min(B,C) >= min(A,D) && max(B,C) <= max(A,D)
             B *<--/          Slope B-C is counted and removed from residue
              / \ /
             /   * C
          \ /
           * A

These steps are fully documented in standards such as
ASTM E 1049 "Standard Practices for Cycle Counting in Fatigue Analysis" [1]
This implementation uses the 4-point algorithm mentioned in [3,4] and the 3-point HCM method proposed in [2] as well as the ASTM E 1049 (2011) standard in [1].
To take the residue into account, you may implement a custom method or use some
predefined functions.

### Features of this package
 1. Modular architecture in two layers:
    a) Module _rainflow.c_ (with _rainflow.h_) holds all necessary functions for rainflow counting and histogram extraction. You may select multiple optional features at compile time:
    `RFC_MINIMAL`: To use core functions for rainflow counting only (for porting to µControllers for example).
    `RFC_TP_SUPPORT`: Turning point storage.
    `RFC_HCM_SUPPORT`: HCM algorithm (Clormann/Seeger).
    `RFC_ASTM_SUPPORT`: ASTM E 1049 (2011) algorithm.
    `RFC_AT_SUPPORT`: User defined amplitude transformation (Haigh diagram).
    `RFC_DH_SUPPORT`:  Damage history storage.
    `RFC_USE_DELEGATES`: Delegates for various core functions to implement user defined behavior.
    `RFC_GLOBAL_EXTREMA`:  Store global data extrema.
    `RFC_DAMAGE_FAST`: Using lookup tables for damage and amplitude transformation.
    `RFC_EXPORT_MEX`: Export a mexFunction() to use the rainflow counting in MATLAB (R).
    `RFC_EXPORT_PY`: Export a Python extension to use the rainflow counting in Python.
    `RFC_UNIT_TEST`: Build an executable for unit testing.
    Using _COAN_ [[http://coan2.sourceforge.net/]] for example, you can tidy up the code from unwanted features.</br>
    (The minimal version of this package is created using _COAN_ with option `RFC_MINMAL`set.)
    b) C++ wrapper _rainflow.hpp_ encapsulates functions from rainflow.h in a namespace and offers a template class _Rainflow_ for object oriented access and inheritance.
    This class also offers container class based turning point storage.
 2. Streamable: You're able to count your data at once, as data packages or sample-wise.
 3. Class width fit to your needs. Dynamically increase class width, when needed. (Needs turning point storage.)
 4. Four point counting method, optionally HCM counting method (Clormann/Seeger).
 5. Woehler curve with up to two slopes, fatigue limit and omission.
 6. Miners rule for damage accumulation (elementary, original, modified and consequent).
 7. In-time damage indicator (Miners' consequent rule).
 8. In-time histograms: rainflow matrix, level crossing and range pair counting.
 9. Turning points with hysteresis filtering. Turning points involved in a closed hysteresis are marked as pairs, with its partial assigned damage. (Compact history)
 10. Look-up tables for damage calculation and amplitude transformation.
 11. Amplitude transformation (Haigh diagram) according to FKM (symmetrical, non-symmetrical or user defined).
 12. Damage history (uncompressed)
 13. Various methods on residual data:
     - According to DIN 45667
     - ASTM method (halfcycle, fullcycle)
     - Second run
     - HCM
 14. Various function pointers to implement user defined behavior.
 15. Conversions supporting RFM->LC, RFM->RP, RFM->Damage and RP->Damage (original, elementar, modifiziert, konsequent).


## Build from sources

Building makefiles from sources root folder in a shell:

    cmake -S. -Bbuild -G "Visual Studio 16 2019"
    cmake --build build --config Release

If you intend to use a special MATLAB(R) installation to build the project
just set an environment variable `Matlab_ROOT_DIR` such as

    export Matlab_ROOT_DIR=/usr/local/MATLAB/R2017b
or use the option to configure the build tree:

    cmake -S. -Bbuild -DMatlab_ROOT_DIR=/usr/local/MATLAB/R2017b -G "Visual Studio 16 2019"
    cmake --build build --config Release

### MATLAB (only)
    cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_UNIT_TEST=0 -G "Visual Studio 16 2019"
    cmake --build build --target rfc_mex --config Release

### Python (only)
    cmake -S. -Bbuild -DRFC_EXPORT_MEX=0 -DRFC_UNIT_TEST=0 -G "Visual Studio 16 2019"
    cmake --build build --target rfcnt --config Release

To build a Python package (wheel) you don't need cmake. Use packages `setuptools`, `wheel` and `build` instead:

    cd src/python
    pip install setuptools build wheel oldest-supported-numpy
    python -mbuild -nw

#### Run examples
    python -m rfcnt.run_examples

Using rfcnt in Google Colaboratory
https://colab.research.google.com

    !pip install --no-build-isolation --no-deps https://github.com/a-ma72/rainflow/releases/download/rfcnt-0.4.7rc0/rfcnt-0.4.7rc0.tar.gz

    import rfcnt
    rfcnt.tests.examples.example_1()

### Unit test (only)
Currently, two options are offered to perform a unit test:
1. Build `rfc_test` and execute in a shell:

        cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_EXPORT_MEX=0 -G "Visual Studio 16 2019"
        cmake --build build --target rfc_test --config Release

    Then invoke the unit test:

        build/test/Release/rfc_test      # Linux
        build\test\Release\rfc_test.exe  # Windows

    or

        ctest -C Release
2. Run Python code (shell):

    python -m rfcnt.run_tests



---
### References:
[1] "Standard Practices for Cycle Counting in Fatigue Analysis."<br>
    ASTM Standard E 1049, 1985 (2011).<br>
    West Conshohocken, PA: ASTM International, 2011.<br>
[2] "Rainflow - HCM / Ein Hysteresisschleifen-Zaehlalgorithmus auf werkstoffmechanischer Grundlage"<br>
    U.H. Clormann, T. Seeger<br>
    1985 TU Darmstadt, Fachgebiet Werkstoffmechanik<br>
[3] "Zaehlverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"<br>
    FVA-Richtlinie, 2010.<br>
    [https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf]<br>
[4] Siemens Product Lifecycle Management Software Inc., 2018.<br>
    [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]<br>
[5] "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"<br>
    G.Marsh;<br>
    International Journal of Fatigue 82 (2016) 757-765,<br>
    [https://doi.org/10.1016/j.ijfatigue.2015.10.007]<br>
[6] "Betriebsfestigkeit - Verfahren und Daten zur Bauteilberechnung"<br>
    Haibach, Erwin; Springer Verlag<br>
[]  "Schaedigungsbasierte Hysteresefilter"; Hack, M, D386 (Diss Univ. Kaiserslautern), Shaker Verlag Aachen, 1998, ISBN 3-8265-3936-2<br>
[]  "Hysteresis and Phase Transition"<br>
    Brokate, M.; Sprekels, J.; Applied Mathematical Sciences 121; Springer, New York, 1996<br>
[]  "Rainflow counting and energy dissipation in elastoplasticity"; Eur. J. Mech. A/Solids 15, pp. 705-737, 1996<br>
    Brokate, M.; Dressler, K.; Krejci, P.<br>
[]  "Multivariate Density Estimation: Theory, Practice and Visualization". New York, Chichester, Wiley & Sons, 1992<br>
    Scott, D.<br>
[]  "Werkstoffmechanik - Bauteile sicher beurteilen undWerkstoffe richtig einsetzen";<br>
     Ralf Buergel, Hans Albert Richard, Andre Riemer; Springer FachmedienWiesbaden 2005, 2014<br>
[]  "Zaehlverfahren und Lastannahme in der Betriebsfestigkeit";<br>
    Michael Koehler, Sven Jenne / Kurt Poetter, Harald Zenner; Springer-Verlag Berlin Heidelberg 2012<br>
