from typing import Optional, Union

from . import ArrayLike, LCMethod, ResidualMethod, RPDamageCalcMethod, SDMethod

def rfc(
    data: ArrayLike,
    class_width: float,
    *,
    class_count: Optional[int] = 100,
    class_offset: Optional[float] = None,
    hysteresis: Optional[float] = None,
    residual_method: Optional[Union[int, ResidualMethod]] = ResidualMethod.REPEATED,
    spread_damage: Optional[Union[int, SDMethod]] = SDMethod.TRANSIENT_23c,
    lc_method: Optional[Union[int, LCMethod]] = LCMethod.SLOPES_UP,
    use_HCM: Optional[Union[int, bool]] = 0,
    use_ASTM: Optional[Union[int, bool]] = 0,
    enforce_margin: Optional[Union[int, bool]] = 0,
    auto_resize: Optional[Union[int, bool]] = 0,
    wl: Optional[dict] = None,
) -> tuple:
    r"""Rainflow counting.

    Parameters
    ----------
    data : ArrayLike
        The input timeseries.
    class_count : Optional[int] = 100
        The number of bins for counting.
        The range of the input series will be spread evenly over all classes (bins).
    class_offset : Optional[float] = 0
        The lower bound of the first bin. If offset and width are set manually,
        take into account that `np.max(data) < class_offset + class_count * class_width`.
        Counting fails otherwise, if `auto_resize` is not set to True.
    class_width : Optional[float] = 1
        The evenly width of each bin.
        If offset and width are set manually, take into account that
        `np.max(data) < class_offset + class_count * class_width`.
        Counting fails otherwise, if `auto_resize` is not set to True.
    hysteresis : Optional[float] = class_width
        The width of the hysteresis filter (also called "Rueckstellbreite" in german).
    residual_method : Optional[Union[int, ResidualMethod]] = ResidualMethod.REPEATED
        How to deal with the residuum (non closed cycles), see ResidualMethod enum for options.
    use_HCM : Optional[Union[int, bool]] = 0
        Whether to use the HCM (Clormann/Seeger) counting method.
    use_ASTM : Optional[Union[int, bool]] = 0
        Whether to use the ASTM counting method.
    enforce_margin : Optional[Union[bool, int]] = True
        Ensuring first and last turning point match to the input timeseries, disregarding hysteresis filtering.
    auto_resize : Optional[Union[bool, int]] = False
        Expand the class range, if value range exceeds while counting. The width of the bins are kept.
    spread_damage : Optional[int] = 8
        How to distribute damage over turning points::

                          * (P4)
             (P2)    *   / (P3c)
                    / \ /
                   /   *    (P3)
             (P1) *

        - 0 = Halfway damage each point (P2,P3)
        - 1 = Damages for linear amplitude ramp over P2 to P3
        - 2 = Damage linear distributed over P2 to P3
        - 3 = Damages for linear amplitude ramp over P2 to P4
        - 4 = Damage linear distributed over P2 to P4
        - 5 = Full damage assigned to P2
        - 6 = Full damage assigned to P3
        - 7 = Damages transient distributed over P2 to P3
        - 8 = Damages transient distributed over P2 to P3c
    lc_method : Optional[int] = 0
        How to count level crossings.
        - 0 = rising slopes only
        - 1 = falling slopes only
        - 2 = rising and falling slopes
    wl: Optional[dict] = dict(sx=1000, nx=1e7, sd=0, nd=np.inf, k=5, k2=k)
        Definition of the SN-curve.

        * `sx` : float
            SN-curve knee between `k` and `k2`.
        * `nx` : float
            Cycle count at `sx`.
        * `sd` : float
            The fatigue strength (SN-curve).
        * `nd` : float
            Cycle count at `sd`.
        * `k` : float
            The slope of the SN-curve for sa >= sx.
        * `k2` : float
            The slope of the SN-curve for sa < sx.
        * `omission` : float
            The omission value. Values where sa < omission are ignored.

    Returns
    -------
    results : dict
        Dictionary containing the following keys:

        - `bkz` : float
            The (pseudo) damage value.

        - `rp` : np.ndarray
            The range pair histogram with shape (cc, 2), where `cc` is the class count.
            The first column contains the range in ascending order (0:cc-1) * cw, and the second column the counts,
            where `cw` is the class width.

        - `lc` : np.ndarray
            The level crossing histogram with shape (cc, 2), where `cc` is the class count.
            The first column contains the class upper levels in ascending order (1:cc) * cw, and the second column the counts,
            where `cw` is the class width.

        - `tp` : np.ndarray
            The turning point information with shape (n, 3), where n is the number of turning points in `data`.
            The first column refers to the turning point as index (1-based) in `data`. The second column contains the
            value (data[index-1]) of the turning point. The third column contains the pseudo damage at this point.
            Note that these pseudo damages contain fractions from other turning points due to `spread_damage` setting.

        - `res_raw` : np.ndarray
            The residuum of the rainflow counting before applying residual methods.

        - `res` : np.ndarray
            The residuum of the rainflow counting after applying residual methods.

        - `rfm` : np.ndarray
            The rainflow matrix of shape (cc, cc), where `cc` is the class count.
            The matrix contains counts for closed cycles over ranges indexed by start-class and stop-class pairs.

        - `dh` : np.ndarray
            The damage history time series, adjacent to the input time series `data`.
            The SN-curve parameters of the "pre-damaged" material after applying the time series (loads).

        - `wl_miner_consistent` : dict
            Dictionary with Miner-consistent SN-curve values:

            * `sx` : float
                SN-curve knee between `k` and `k2`.
            * `nx` : float
                Cycle count at `sx`.
            * `sd` : float
                The fatigue strength (SN-curve).
            * `nd` : float
                Cycle count at `sd`.
            * `k` : float
                The slope of the SN-curve for sa >= sx.
            * `k2` : float
                The slope of the SN-curve for sa < sx.
            * `q` : float
                Degradation parameter at `sx`, `nx`.
            * `q2` : float
                Degradation parameter at `sd`, `nd`.
            * `omission` : float
                The omission value. Values where sa < omission are ignored.
            * `D` : float
                (Pseudo) damage value according to the SN-curve and
                Miner's consistent rule.

    """

def damage_from_rp(
        Sa: ArrayLike,
        counts: ArrayLike,
        *,
        wl: Optional[dict] = None,
        method: Optional[RPDamageCalcMethod] = 0,
) -> float:
    r"""Calculate damage value from range pair histogram.

    Parameters
    ----------
    Sa : ArrayLike
        Amplitude vector.
    counts : ArrayLike
        Cycle counts vector.
    wl: Optional[dict] = dict(sx=1000, nx=1e7, sd=0, nd=np.inf, k=5, k2=k)
        Definition of the SN-curve.

        * `sx` : float
            SN-curve knee between `k` and `k2`.
        * `nx` : float
            Cycle count at `sx`.
        * `sd` : float
            The fatigue strength (SN-curve).
        * `nd` : float
            Cycle count at `sd`.
        * `k` : float
            The slope of the SN-curve for sa >= sx.
        * `k2` : float
            The slope of the SN-curve for sa < sx.
        * `omission` : float
            The omission value. Values where sa < omission are ignored.
    method : Optional[int] = 0
        Method for damage calculation.

        - 0 = Default (by given `wl`)
        - 1 = Miner elementar
        - 2 = Miner modified
        - 3 = Miner consistent

    Returns
    -------
    damage : float
        (Pseudo) damage value.

    """
