import sys
import unittest
from io import StringIO
from pprint import pprint

from .tests import test_rfcnt


def run() -> unittest.result.TestResult:
    stream = StringIO()
    runner = unittest.TextTestRunner(stream=stream)
    test_result = runner.run(unittest.makeSuite(test_rfcnt.TestRainflowCounting))
    print('Tests run ', test_result.testsRun)
    print('Errors ', test_result.errors)
    pprint(test_result.failures)
    stream.seek(0)
    print('Test output\n', stream.read())
    return test_result


if __name__ == '__main__':
    result = run()
    if result.wasSuccessful():
        sys.exit(0)
    else:
        sys.exit(1)
