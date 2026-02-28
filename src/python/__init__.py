"""Rainflow Counting Package.

This package provides tools and enumerations for rainflow cycle counting and fatigue analysis,
including methods for handling residuals, spreading damage, and level crossing counting.
"""

from __future__ import annotations

from collections import namedtuple
from enum import IntEnum

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
    """An enumeration which slopes encounter level crossing counting.

    Attributes
    ----------
    SLOPES_UP : int
        Count on rising slopes only (default).
    SLOPES_DOWN : int
        Count on falling slopes only.
    SLOPES_ALL : int
        Count on rising AND falling slopes.

    """
    SLOPES_UP = 0           # Count on rising slopes only (default).
    SLOPES_DOWN = 1         # Count on falling slopes only.
    SLOPES_ALL = 3          # Count on rising AND falling slopes.


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


# Import python extension
from . import rfcnt  # noqa: E402 I001 F401

# For backward compatibility, supporting both rfcnt.rfc() and rfcnt.rfcnt.rfc()
from .rfcnt import rfc, damage_from_rp  # noqa: E402 F401

# from . import tests, utils  # noqa: F402
del annotations, NumpyVersion, namedtuple
