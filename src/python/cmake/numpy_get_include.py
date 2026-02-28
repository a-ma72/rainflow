"""Module to generate a header include path for numpy's arrayobject.h."""

import logging
from pathlib import Path

import numpy as np


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

numpy_include_dir = np.get_include()
include_file = (Path(numpy_include_dir) / "numpy" / "arrayobject.h").as_posix()
with Path("numpy_arrayobject.h").open("wt") as f:
    f.write(f'#include "{include_file}"')

logger.info(Path(numpy_include_dir).as_posix())