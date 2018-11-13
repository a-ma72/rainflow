# Rainflow Counting Algorithm (4-point-method), C99 compliant
  
"Rainflow Counting" consists of four main steps:
  1. Hysteresis Filtering
  2. Peak-Valley Filtering
  3. Discretization
  4. Four Point Counting Method:
```
                     * D
                    / \       Closed, if min(B,C) >= min(A,D) && max(B,C) <= max(A,D)
             B *<--/          Slope B-C is counted and removed from residue
              / \ /
             /   * C
          \ /
           * A
```
These steps are fully documented in standards such as 
ASTM E1049 "Standard Practices for Cycle Counting in Fatigue Analysis" [1].
This implementation uses the 4-point algorithm mentioned in [3,4] and the 3-point HCM method proposed in [2].
To take the residue into account, you may implement a custom method or use some
predefined functions.
 
---
### References:
[1] ASTM Standard E 1049, 1985 (2011).   
    "Standard Practices for Cycle Counting in Fatigue Analysis."  
    West Conshohocken, PA: ASTM International, 2011.  
[2] Rainflow - HCM  
    "Ein Hysteresisschleifen-Zaehlalgorithmus auf werkstoffmechanischer Grundlage"  
    U.H. Clormann, T. Seeger  
    1985 TU Darmstadt, Fachgebiet Werkstoffmechanik  
[3] FVA-Richtlinie, 2010.  
    "Zaehlverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"  
    [https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf]  
[4] Siemens Product Lifecycle Management Software Inc., 2018.   
    [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]  
[5] G.Marsh on: "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"  
    International Journal of Fatigue 82 (2016) 757-765,  
    [https://doi.org/10.1016/j.ijfatigue.2015.10.007]  
[]  Hack, M: Schaedigungsbasierte Hysteresefilter; D386 (Diss Univ. Kaiserslautern), Shaker Verlag Aachen, 1998, ISBN 3-8265-3936-2  
[]  Brokate, M; Sprekels, J, Hysteresis and Phase Transition, Applied Mathematical Sciences 121, Springer,  New York, 1996  
[]  Brokate, M; Dressler, K; Krejci, P: Rainflow counting and energy dissipation in elastoplasticity, Eur. J. Mech. A/Solids 15, pp. 705-737, 1996  
[]  Scott, D.: Multivariate Density Estimation: Theory, Practice and Visualization. New York, Chichester, Wiley & Sons, 1992  
