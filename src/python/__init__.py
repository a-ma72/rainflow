from __future__ import annotations

from collections import namedtuple
from enum import IntEnum
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


def _npy_version() -> "NumpyVersion":
    """
    Query the current version of NumPy.

    This function retrieves the current version of NumPy and returns it as a
    NumpyVersion object.

    Returns
    -------
    NumpyVersion
        The current NumPy version.

    Notes
    -----
    For more information, see the official documentation:
    https://numpy.org/doc/stable/reference/generated/numpy.lib.NumpyVersion.html

    Examples
    --------
    >>> _npy_version()
    NumpyVersion('1.21.0')
    """
    from numpy import __version__ as version  # Import NumPy's version string
    from numpy.lib import NumpyVersion  # Import the NumpyVersion class

    return NumpyVersion(version)  # Return the version as a NumpyVersion object


def _npy_get_capi(npy_version: "NumpyVersion") -> int:
    """
    Return the NumPy C API VERSION based on the given NumPy version.

    Parameters
    ----------
    npy_version : NumpyVersion
        The NumPy version to check against.

    Returns
    -------
    int
        The corresponding NumPy C API version.

    Notes
    -----
    You can get the C API Version of your current NumPy installation by inspecting
    the `NPY_x_y_API_VERSION` in "numpy/core/include/numpy/numpyconfig.h".
    NumPy's "C API Version" does not increment with every release.

    The mapping from NumPy version to C API version is as follows:
    - NumPy < 1.20: returns 0x0d (fallback to NumPy 1.19)
    - NumPy < 1.22: returns 0x0e (fallback to NumPy 1.21)
    - NumPy < 1.23: returns 0x0f (fallback to NumPy 1.22)
    - NumPy < 2.0:  returns 0x10 (fallback to NumPy 1.24)
    - NumPy >= 2.0: returns 0x12 (NumPy 2.0 and above)

    Examples
    --------
    >>> from numpy.lib import NumpyVersion
    >>> _npy_get_capi(NumpyVersion('1.21.0'))
    14
    >>> _npy_get_capi(NumpyVersion('2.0.0'))
    18

    """
    if npy_version < '1.20':
        return 0x0d  # Fallback to NumPy 1.19
    elif npy_version < '1.22':
        return 0x0e  # Fallback to NumPy 1.21
    elif npy_version < '1.23':
        return 0x0f  # Fallback to NumPy 1.22
    elif npy_version < '2.0':
        return 0x10  # Fallback to NumPy 1.24
    else:
        return 0x12  # NumPy 2.0 and above


def _load_spec_extension_module():
    """
    Select and load a suitable extension build based on the NumPy version and Python version.

    This function dynamically selects an appropriate extension module for the current
    environment, ensuring compatibility with the specific NumPy version and Python version
    in use.

    Raises
    ------
    ImportError
        If no suitable extension build is found.

    Notes
    -----
    - The function requires Python 3.8 or higher.
    - The extension module is expected to be located in a directory named `_ext` within the package.

    Examples
    --------
    >>> _load_spec_extension_module()
    """
    import os
    import sys
    from importlib.util import module_from_spec, spec_from_file_location

    EXT_DIR = "_ext"  # Directory containing the extension modules

    # Ensure Python version is 3.8 or higher
    assert sys.version_info >= (3, 8), "Only Python >= 3.8 supported."

    prebuilds = [file for file in os.listdir(os.path.join(__path__[0], EXT_DIR)) if file.endswith(".pyd")]

    for prebuild in prebuilds:
        spec = spec_from_file_location(".rfcnt", os.path.join(__path__[0], EXT_DIR, prebuild))
        if spec:
            try:
                module = module_from_spec(spec)
                spec.loader.exec_module(module)
                sys.modules[__name__ + spec.name] = module
                print(f"Using prebuild `{prebuild}`")
                break
            except:
                pass
    else:
        npy_version = _npy_version()  # Get the current NumPy version
        raise ImportError(f"No suitable rfcnt build found for NumPy {npy_version.vstring}!")


if _npy_version() >= "1.20.0":
    from numpy.typing import ArrayLike
else:
    from typing import Any as ArrayLike


# Try to import the extension module, compiled from sources.
# When this module is installed from a wheel, there are multiple
# versions depending on the version of NumPy installed.
# In the latter case, the suitable version will be loaded and
# finally imported.
try:
    from . import rfcnt
except (ImportError, RuntimeError):
    import os
    if os.name == "nt":
        # If the operating system is MS Windows, try to load a suitable prebuilt module.
        _load_spec_extension_module()
        from . import rfcnt
        print("Success: Loaded another suitable prebuilt instead.")
    else:
        # For non-Windows systems, re-raise the exception if loading fails.
        raise

# For backward compatibility, supporting both rfcnt.rfc() and rfcnt.rfcnt.rfc()
rfc = rfcnt.rfc

from . import tests, utils  # noqa F402
