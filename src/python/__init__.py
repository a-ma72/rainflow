from __future__ import annotations

import json
import os
import warnings
from collections import namedtuple
from enum import IntEnum
from numpy import __version__ as _npy_version
from numpy.lib import NumpyVersion
from .version import __version__


ClassParams = namedtuple("ClassParams", "class_count, class_offset, class_width")


class ResidualMethod(IntEnum):
    """
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
    """
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
    """
    An enumeration which slopes encounter level crossing counting.

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
    """
    A method enumeration how `damage_from_rp()` calculates the damage value.

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


def _get_spec_extension_prebuild():
    if os.name == "nt":
        import sys
        from importlib.util import module_from_spec, spec_from_file_location

        EXT_DIR = "_ext"  # Directory containing the extension modules

        # Ensure Python version is 3.8 or higher
        if sys.version_info < (3, 8):
            warnings.warn("Prebuilds are supported for Python >= 3.8 only.")
            return

        # Determine the directory of the current script
        package_directory = os.path.dirname(__file__)

        # Try to locate suitable prebuilt modules.
        prebuilds = None
        if os.path.exists(os.path.join(package_directory, "_ext")):
            # Construct the path to the text file
            prebuilds_json_path = os.path.join(package_directory, "prebuilds.json")
            try:
                with open(prebuilds_json_path, "rt") as f:
                    prebuilds = json.load(f)
            except Exception:
                pass

        npy_version = NumpyVersion(_npy_version)
        candidates = []
        if prebuilds:
            for version_str in prebuilds:
                version = NumpyVersion(version_str)
                if version.major == npy_version.major and version <= npy_version:
                    candidates.append(prebuilds[version_str]["target"])

        if candidates and "root" not in candidates:
            prebuild = candidates[0]
            # TODO: Support other Python tags too?
            files = [file for file in os.listdir(os.path.join(package_directory, EXT_DIR)) if file.startswith(prebuild)]

            for file in files:
                spec = spec_from_file_location(".rfcnt", os.path.join(package_directory, EXT_DIR, file))
                if spec:
                    try:
                        module = module_from_spec(spec)
                        spec.loader.exec_module(module)
                        sys.modules[__name__ + spec.name] = module
                    except:
                        raise ImportError(f"No suitable rfcnt build found for NumPy {_npy_version}!")
                return prebuild


if NumpyVersion(_npy_version) >= "1.20.0":
    from numpy.typing import ArrayLike
else:
    from typing import Any as ArrayLike


# Try to import the extension module, compiled from sources.
# When this module is installed from a wheel, there are multiple
# versions depending on the version of NumPy installed.
# In the latter case, the suitable version will be loaded and
# finally imported.
_prebuild = _get_spec_extension_prebuild()
if _prebuild:
    print(f"Using prebuild `{_prebuild}`")

# Import python extension
from . import rfcnt  # noqa 402

# For backward compatibility, supporting both rfcnt.rfc() and rfcnt.rfcnt.rfc()
from .rfcnt import rfc, damage_from_rp  # noqa 402

# from . import tests, utils  # noqa F402
del _get_spec_extension_prebuild, annotations, NumpyVersion, namedtuple, os, json, warnings
