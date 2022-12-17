import os
import sys
from enum import Enum
from importlib.util import module_from_spec, spec_from_file_location
from typing import Union, Optional

from numpy.typing import ArrayLike
from numpy import asarray_chkfinite, double as np_double, pad as np_pad


def _load_spec_extension_module():
    """Select and load a suitable extension build"""
    from numpy import __version__ as np__version__
    from numpy.lib import NumpyVersion

    EXT_DIR = "_ext"

    # see `numpy/numpy/core/include/numpy/numpyconfig.h`
    npy_version = NumpyVersion(np__version__)
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
        raise ImportError("No suitable build found for numpy %s!", np__version__)

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


class ResidualMethod(Enum):
    NONE = 0
    _IGNORE = 1
    _NO_FINALIZE = 2
    DISCARD = 3
    HALFCYCLES = 4
    FULLCYCLES = 5
    CLORMANN_SEEGER = 6
    REPEATED = 7
    DIN45667 = 8

class SDMethod(Enum):
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

class LCMethod(Enum):
    SLOPES_UP = 0
    SLOPES_DOWN = 1
    SLOPES_ALL = 3


def rfc(data: ArrayLike,
        class_count: Optional[int] = 100,
        class_width: Optional[float] = None,
        class_offset: Optional[float] = None,
        hysteresis: Optional[float] = None,
        residual_method: Optional[Union[int, ResidualMethod]] = ResidualMethod.REPEATED,
        spread_damage: Optional[Union[int, SDMethod]] = SDMethod.TRANSIENT_23c,
        lc_method: Optional[Union[int, LCMethod]] = LCMethod.SLOPES_UP,
        use_HCM: Optional[Union[int, bool]] = 0,
        use_ASTM: Optional[Union[int, bool]] = 0,
        enforce_margin: Optional[Union[int, bool]] = 0,
        auto_resize: Optional[Union[int, bool]] = 0,
        wl: Optional[dict] = None) -> tuple:
    """Wrapper for .rfcnt.rfc()
    """
    data = asarray_chkfinite(data, dtype=np_double, order="C").flatten()
    if class_width is None:
        class_width = data.ptp() / (class_count - 1)
    if class_offset is None:
        class_offset = data.min() - class_width / 2
    if hysteresis is None:
        hysteresis = class_width
    if wl is None:
        wl = dict(sd=1e3, nd=1e7, k=5, k2=5)
    if isinstance(residual_method, ResidualMethod):
        residual_method = residual_method.value
    if isinstance(spread_damage, SDMethod):
        spread_damage = spread_damage.value
    if isinstance(lc_method, LCMethod):
        lc_method = lc_method.value
    res = rfcnt.rfc(
        data,
        class_count=class_count,
        class_offset=class_offset,
        class_width=class_width,
        hysteresis=hysteresis,
        use_HCM=int(use_HCM),
        use_ASTM=int(use_ASTM),
        enforce_margin=int(enforce_margin),
        auto_resize=int(auto_resize),
        spread_damage=spread_damage,
        residual_method=residual_method,
        lc_method=lc_method,
        wl=wl)
    if res["dh"].size < data.size:
        res["dh"] = np_pad(res["dh"], (0, data.size - res["dh"].size))
    return res


from . import tests, utils
