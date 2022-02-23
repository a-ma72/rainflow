from rfcnt.tests import test_rfcnt
import unittest
from io import StringIO
from pprint import pprint


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


"""
Jupyter Notebook:
!pip install ./rfcnt-0.3.1.tar.gz
!pip install --upgrade matplotlib
!python -m rfcnt.run_tests
"""
