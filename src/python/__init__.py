"""Rainflow Counting Package.

This package provides tools and enumerations for rainflow cycle counting and fatigue
analysis, including methods for handling residuals, spreading damage, and level
crossing counting.
"""

from __future__ import annotations

import os
from collections import namedtuple
from enum import IntEnum
from pathlib import Path
from typing import Union

from numpy import __version__ as _npy_version
from numpy.lib import NumpyVersion

from .version import __version__

ClassParams = namedtuple("ClassParams", "class_count, class_offset, class_width")


class ResidualMethod(IntEnum):
    """Enum for residual methods in rainflow counting.

    An enumeration representing various methods for handling residuals in data analysis.

    Attributes
    ----------
    NONE : int
        No residual method applied.
    _IGNORE : int
        Ignore residuals in the computation.
    _NO_FINALIZE : int
        Do not finalize the computation with residuals.
    DISCARD : int
        Discard residuals completely.
    HALFCYCLES : int
        Apply the half cycles method for residuals.
    FULLCYCLES : int
        Apply the full cycles method for residuals.
    CLORMANN_SEEGER : int
        Use the Clormann-Seeger method for residuals.
    REPEATED : int
        Use repeated application of residuals method.
    DIN45667 : int
        Apply the DIN 45667 standard method for residuals.

    """

    NONE = 0             # No residual method applied.
    _IGNORE = 1          # Ignore residuals in the computation.
    _NO_FINALIZE = 2     # Do not finalize the computation with residuals.
    DISCARD = 3          # Discard residuals completely.
    HALFCYCLES = 4       # Apply the half cycles method for residuals.
    FULLCYCLES = 5       # Apply the full cycles method for residuals.
    CLORMANN_SEEGER = 6  # Use the Clormann-Seeger method for residuals.
    REPEATED = 7         # Use repeated application of residuals method.
    DIN45667 = 8         # Apply the DIN 45667 standard method for residuals.


class SDMethod(IntEnum):
    """An enumeration for methods of spreading damage.

    An enumeration representing various methods spreading damage increments over time history.

    Attributes
    ----------
    NONE : int
        No spread damage calculation.
    HALF_23 : int
        Equally split damage between P2 and P3.
    RAMP_AMPLITUDE_23 : int
        Spread damage according to amplitude over points between P2 and P3.
    RAMP_DAMAGE_23 : int
        Spread damage evenly over points between P2 and P3.
    RAMP_AMPLITUDE_24 : int
        Spread damage exponentially according to amplitude impact over points between P2 and P4.
    RAMP_DAMAGE_24 : int
        Spread damage evenly over points between P2 and P4.
    FULL_P2 : int
        Assign damage to P2.
    FULL_P3 : int
        Assign damage to P3.
    TRANSIENT_23 : int
        Spread damage transient according to amplitude over points between P2 and P3.
    TRANSIENT_23c : int
        Spread damage transient according to amplitude over points between P2 and P4 only until cycle is closed.

    """

    NONE = -1               # No spread damage calculation.
    HALF_23 = 0             # Equally split damage between P2 and P3.
    RAMP_AMPLITUDE_23 = 1   # Spread damage according to amplitude over points between P2 and P3.
    RAMP_DAMAGE_23 = 2      # Spread damage evenly over points between P2 and P3.
    RAMP_AMPLITUDE_24 = 3   # Spread damage exponentially according to amplitude impact over points between P2 and P4.
    RAMP_DAMAGE_24 = 4      # Spread damage evenly over points between P2 and P4.
    FULL_P2 = 5             # Assign damage to P2.
    FULL_P3 = 6             # Assign damage to P3.
    TRANSIENT_23 = 7        # Spread damage transient according to amplitude over points between P2 and P3.
    TRANSIENT_23c = 8       # Spread damage transient according to amplitude over points between P2 and P4 only until cycle is closed.


class LCMethod(IntEnum):
    """Which slopes contribute to level-crossing counts.

    Simple enumeration matching ``rfc_lc_count_method`` (not a bit mask).
    ``0 | 1`` is DOWN, not ALL.

    DIN 45667 uses a **static** crossing direction for every class bound
    (``SLOPES_UP``, ``SLOPES_DOWN``, or ``SLOPES_ALL``). FVA Merkblatt uses a
    **sign-dependent** direction from a zero-load baseline (``FVA``).

    Attributes
    ----------
    SLOPES_UP : int
        DIN 45667: count on rising slopes only (value 0).
    SLOPES_DOWN : int
        DIN 45667: count on falling slopes only (value 1).
    SLOPES_ALL : int
        DIN 45667: count on rising AND falling slopes (value 2). Default for
        ``rfc()`` and ``RFC``, matching C ``RFC_FLAGS_COUNT_LC``.
    FVA : int
        FVA Merkblatt: count both slopes internally, then convert on read so
        positive-going crossings apply for class upper bounds ``u >= 0`` and
        negative-going for ``u < 0``. Residue methods do not change this
        histogram.
    DIN45667 : int
        Compatibility alias of :attr:`FVA` (historical misnomer).

    """

    SLOPES_UP = 0           # DIN 45667: rising slopes only (static global direction).
    SLOPES_DOWN = 1         # DIN 45667: falling slopes only (static global direction).
    SLOPES_ALL = 2          # DIN 45667: rising AND falling slopes (Python/C default).
    FVA = 3                 # FVA Merkblatt: sign-dependent (UP for u>=0, DOWN for u<0).
    DIN45667 = 3            # Compatibility alias of FVA (historical misnomer).


class RPDamageCalcMethod(IntEnum):
    """A method enumeration how `damage_from_rp()` calculates the damage value.

    Attributes
    ----------
    DEFAULT : int
        Use SN curve params as they are set.
    MINER_ELEMENTAR : int
        Use SN curve type "Miner elementar".
        (Slope `k2` ignored.)
    MINER_MODIFIED : int
        Use SN curve type "Miner modified".
        (Takes slope `k2` into account.)
    MINER_CONSISTENT : int
        Accumulate according to "Miner consistent".

    """

    DEFAULT = 0             # Use SN curve params as they are set.
    MINER_ELEMENTAR = 1     # Use SN curve type "Miner elementary".
    MINER_MODIFIED = 2      # Use SN curve type "Miner modified".
    MINER_CONSISTENT = 3    # Accumulate according to "consistent Miner's rule".


if NumpyVersion(_npy_version) >= "1.20.0":
    from numpy.typing import ArrayLike
else:
    from typing import Any as ArrayLike


def _add_vendored_dll_dir() -> None:
    """Add delvewheel ``rfcnt.libs`` to the Windows DLL search path.

    Wheels ship ``msvcp140-*.dll`` next to the package. Python 3.8+ does not
    search that folder unless ``os.add_dll_directory`` is called first.
    """
    if os.name != "nt":
        return
    libs = Path(__file__).resolve().parent.parent / "rfcnt.libs"
    if not libs.is_dir():
        return
    add_dir = getattr(os, "add_dll_directory", None)
    if add_dir is not None:
        add_dir(str(libs))


_add_vendored_dll_dir()

# Import python extension
from . import rfcnt  # noqa: E402 I001 F401

# For backward compatibility, supporting both rfcnt.rfc() / rfcnt.RFC()
# and rfcnt.rfcnt.rfc() / rfcnt.rfcnt.RFC()
from .rfcnt import rfc, damage_from_rp, at_transform  # noqa: E402 F401
from .rfcnt import RFC as _RFC_C  # noqa: E402 F401


class RFC:
    """Stateful rainflow counter.

    Construct with class parameters, then call :meth:`feed` one or more times
    and :meth:`finalize` when the series is complete. The C implementation
    (``rfcnt.rfcnt.RFC``) owns the counting context.

    Damage history (``spread_damage`` other than :data:`SDMethod.NONE`) is not
    supported here; use :func:`rfc` for one-shot counting with a damage history.
    """

    def __init__(self, class_width: float, **kwargs):
        self._impl = _RFC_C(class_width, **kwargs)

    def feed(self, data: ArrayLike) -> RFC:
        """Append samples and continue counting. Returns self for chaining."""
        self._impl.feed(data)
        return self

    def finalize(self, residual_method: Union[int, ResidualMethod] = ResidualMethod.REPEATED) -> None:
        """Close open cycles with the given residual method."""
        return self._impl.finalize(residual_method)

    def damage_as(self, residual_method: Union[int, ResidualMethod] = ResidualMethod.REPEATED) -> float:
        """Return damage as if :meth:`finalize` had been called. Does not change this object."""
        return self._impl.damage_as(residual_method)

    def rp_as(self, residual_method: Union[int, ResidualMethod] = ResidualMethod.REPEATED):
        """Return range pairs as if :meth:`finalize` had been called. Does not change this object."""
        return self._impl.rp_as(residual_method)

    def lc_as(self, residual_method: Union[int, ResidualMethod] = ResidualMethod.REPEATED):
        """Return the level-crossing histogram as if :meth:`finalize` had been called.

        DIN 45667 methods (``SLOPES_UP`` / ``SLOPES_DOWN`` / ``SLOPES_ALL``)
        keep the static global direction. ``FVA`` converts to sign-dependent
        counts (UP for ``u >= 0``, DOWN for ``u < 0``). Does not change this
        object.
        """
        return self._impl.lc_as(residual_method)

    def rfm_as(self, residual_method: Union[int, ResidualMethod] = ResidualMethod.REPEATED):
        """Return the rainflow matrix as if :meth:`finalize` had been called. Does not change this object."""
        return self._impl.rfm_as(residual_method)

    def at_init(self, M: float, **kwargs) -> None:
        """Initialize amplitude transformation (Haigh / FKM) on this counter.

        Must be called after construction and before :meth:`feed`. See
        :func:`at_transform` for parameter names (`R_rig`, `Sm_rig`,
        `R_pinned`, `Sa_ref`, `Sm_ref`, `symmetric`).
        """
        return self._impl.at_init(M, **kwargs)

    def at_transform(self, Sa: ArrayLike, Sm: ArrayLike):
        """Apply the Haigh transformation configured by :meth:`at_init`.

        `Sa` and `Sm` may be scalars or arrays of equal length. Returns
        transformed amplitudes with the same shape as `Sa`.
        """
        return self._impl.at_transform(Sa, Sm)

    def close(self):
        """Release the rainflow context."""
        return self._impl.close()

    def __enter__(self):
        self._impl.__enter__()
        return self

    def __exit__(self, exc_type, exc, tb):
        return self._impl.__exit__(exc_type, exc, tb)

    @property
    def state(self) -> int:
        return self._impl.state

    @property
    def error(self) -> int:
        return self._impl.error

    @property
    def damage(self) -> float:
        return self._impl.damage

    @property
    def residue(self):
        return self._impl.residue

    @property
    def rp(self):
        return self._impl.rp

    @property
    def lc(self):
        return self._impl.lc

    @property
    def rfm(self):
        return self._impl.rfm

    @property
    def tp(self):
        """Turning points as an ``(n, 4)`` array: pos, value, damage, adj_pos."""
        return self._impl.tp

    @property
    def res_raw(self):
        """Open residue before the residual method (4-point strip; isolated read-only snapshot after finalize)."""
        return self._impl.res_raw

    @property
    def wl_miner_consistent(self):
        """Live Miner-consistent (impaired) Wöhler dict; matches ``rfc()["wl_miner_consistent"]`` after finalize()."""
        return self._impl.wl_miner_consistent


# from . import tests, utils  # noqa: F402
del annotations, NumpyVersion, namedtuple
