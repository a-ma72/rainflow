import unittest
from io import StringIO
from pprint import pprint

from .tests import test_rfcnt


def run():
    stream = StringIO()
    runner = unittest.TextTestRunner(stream=stream)
    result = runner.run(unittest.makeSuite(test_rfcnt.TestRainflowCounting))
    print('Tests run ', result.testsRun)
    print('Errors ', result.errors)
    pprint(result.failures)
    stream.seek(0)
    print('Test output\n', stream.read())


if __name__ == '__main__':
    run()
