from enum import IntEnum


class ResidualMethod(IntEnum):
    NONE = 0
    _IGNORE = 1
    _NO_FINALIZE = 2
    DISCARD = 3
    HALFCYCLES = 4
    FULLCYCLES = 5
    CLORMANN_SEEGER = 6
    REPEATED = 7
    DIN45667 = 8


class SDMethod(IntEnum):
    NONE = -1
    HALF_23 = 0
    RAMP_AMPLITUDE_23 = 1
    RAMP_DAMAGE_23 = 2
    RAMP_AMPLITUDE_24 = 3
    RAMP_DAMAGE_24 = 4
    FULL_P2 = 5
    FULL_P3 = 6
    TRANSIENT_23 = 7
    TRANSIENT_23c = 8


class LCMethod(IntEnum):
    SLOPES_UP = 0
    SLOPES_DOWN = 1
    SLOPES_ALL = 3


def _npy_version():
    """Query numpy version.
    (See <https://numpy.org/doc/stable/reference/generated/numpy.lib.NumpyVersion.html>)
    """
    from numpy import __version__ as version
    from numpy.lib import NumpyVersion

    return NumpyVersion(version)


def _load_spec_extension_module():
    """Select and load a suitable extension build"""
    import os, sys
    from importlib.util import module_from_spec, spec_from_file_location
    EXT_DIR = "_ext"

    assert sys.version_info >= (3, 9)

    npy_version = _npy_version()
    if npy_version < '1.20.0':
        # API version 0x8
        infix = "npy_1_19_5"
    elif npy_version < '1.22.0':
        # API version 0xe
        infix = "npy_1_21_6"
    elif npy_version < '1.23.0':
        # API version 0xf
        infix = "npy_1_22_4"
    else:
        # API version 0x10
        infix = "npy_1_23_5"

    modulename = [file for file in os.listdir(os.path.join(__path__[0], EXT_DIR)) if infix in file]

    if len(modulename) == 1:
        spec = spec_from_file_location(".rfcnt", os.path.join(__path__[0], EXT_DIR, modulename[0]))
        module = module_from_spec(spec)
        sys.modules[__name__ + spec.name] = module
        spec.loader.exec_module(module)
        del spec, module, npy_version, infix
    else:
        raise ImportError("No suitable build found for numpy %s!",
                          npy_version.vstring)


if _npy_version() >= "1.20.0":
    from numpy.typing import ArrayLike
else:
    from typing import Any as ArrayLike


# Try to import the extension module, compiled from sources.
# When this module is installed from a wheel, there are multiple 
# versions depending on the version of numpy installed.
# In latter case, the suitable version will be loaded and
# finally imported.
try:
    from . import rfcnt
except ImportError:
    _load_spec_extension_module()
    from . import rfcnt
# For backward compatibility supporting rfcnt.rfc() and rfcnt.rfcnt.rfc()
rfc = rfcnt.rfc
from . import tests, utils  # noqa F402
