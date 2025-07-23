"""Run examples.

Module to run example functions from the examples module for demonstration and testing purposes.
"""

import logging

from .tests import examples

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run example 1 from the examples module."""
    logger.info("Running example 1...")
    examples.example_1()


if __name__ == "__main__":
    main()


"""
Jupyter Notebook:
!pip install ./rfcnt-0.4.7.tar.gz
!pip install -U matplotlib
!python -m rfcnt.run_examples
"""
