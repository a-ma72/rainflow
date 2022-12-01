import os
import sys
from importlib.util import spec_from_file_location, module_from_spec
from . import utils
from . import tests


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
