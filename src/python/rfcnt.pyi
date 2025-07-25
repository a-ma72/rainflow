from typing import Optional, Union

from . import ArrayLike, LCMethod, ResidualMethod, SDMethod, RPDamageCalcMethod


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
        wl: Optional[dict] = None
) -> tuple: ...


def damage_from_rp(
        Sa: ArrayLike,
        counts: ArrayLike,
        *,
        wl: Optional[dict] = None,
        method: Optional[RPDamageCalcMethod] = 0
) -> float: ...
