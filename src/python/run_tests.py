"""Module to run the TestRainflowCounting test suite and log the results.

This script sets up logging, runs the unit tests for the rainflow counting implementation,
and exits with an appropriate status code based on the test results.
"""

import logging
import sys
import unittest
from io import StringIO

from .tests import test_rfcnt

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run() -> unittest.result.TestResult:
    """Run the TestRainflowCounting test suite and prints the results.

    Returns
    -------
    unittest.result.TestResult
        The result object containing information about the test run.

    """
    stream = StringIO()
    runner = unittest.TextTestRunner(stream=stream)
    test_result = runner.run(unittest.makeSuite(test_rfcnt.TestRainflowCounting))
    logger.info("Tests run %d", test_result.testsRun)
    logger.info("Errors %s", test_result.errors)
    logger.info("Failures: %s", test_result.failures)
    stream.seek(0)
    logger.info("Test output\n%s", stream.read())
    return test_result


if __name__ == "__main__":
    result = run()
    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)



"""
Jupyter Notebook:
!pip install ./rfcnt-0.4.7.tar.gz
!pip install -U matplotlib
!python -m rfcnt.run_tests
"""
