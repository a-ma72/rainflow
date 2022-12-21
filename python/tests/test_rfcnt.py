import os
import unittest

import numpy as np

from .. import rfc, ResidualMethod, SDMethod


class TestRainflowCounting(unittest.TestCase):

    def get_script_path(self):
        return os.path.dirname(__file__)

    def class_param(self, data, class_count):
        assert class_count > 1

        if len(data) == 0:
            class_width  = 1  # noqa E221
            class_offset = 0  # noqa E221
        else:
            class_width  = data.ptp() / (class_count - 1)  # noqa E221
            class_width  = np.ceil(class_width * 100) / 100  # noqa E221
            class_offset = np.floor((data.min() - class_width / 2) * 1000) / 1000  # noqa E221, E501

        return class_width, class_offset

    def test_empty_series(self):
        class_count       =  100  # noqa E221
        x                 =  np.array([])  # noqa E221
        x_max             =  1  # noqa E221
        x_min             = -1  # noqa E221
        class_width, \
         class_offset     =  self.class_param(x, class_count)  # noqa E221
        hysteresis        =  class_width  # noqa E221
        enforce_margin    =  False                # First and last data point may be excluded in tp  # noqa E221
        use_HCM           =  False                # Use 4 point method, not HCM  # noqa E221
        use_ASTM          =  False                # Use 4 point method, not ASTM  # noqa E221
        residual_method   =  ResidualMethod.NONE  # No processing on residue  # noqa E221
        spread_damage     =  SDMethod.NONE        # No damage spreading  # noqa E221

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        self.assertEqual(res["rfm"].sum(), 0)
        self.assertEqual(len(res["res"]), 0)

    def test_single_cycle_up(self):
        class_count       =  4  # noqa E221
        x                 =  np.array([1, 3, 2, 4])  # noqa E221
        x_max             =  4  # noqa E221
        x_min             =  1  # noqa E221
        class_width, \
         class_offset     =  self.class_param(x, class_count)  # noqa E221
        hysteresis        =  class_width * 0.99  # noqa E221
        enforce_margin    =  False                # First and last data point may be excluded in tp  # noqa E221
        use_HCM           =  False                # Use 4 point method, not HCM  # noqa E221
        use_ASTM          =  False                # Use 4 point method, not ASTM  # noqa E221
        residual_method   =  ResidualMethod.NONE  # No processing on residue  # noqa E221
        spread_damage     =  SDMethod.NONE        # No damage spreading  # noqa E221

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        self.assertEqual(res["rfm"].sum(), 1)
        self.assertEqual(res["rfm"][3 - 1, 2 - 1], 1)
        self.assertTrue((res["res"].flatten() == [1, 4]).all())

    def test_one_cycle_down(self):
        class_count       =  4  # noqa E221
        x                 =  np.array([4, 2, 3, 1])  # noqa E221
        x_max             =  4  # noqa E221
        x_min             =  1  # noqa E221
        class_width, \
         class_offset     =  self.class_param(x, class_count)  # noqa E221
        hysteresis        =  class_width * 0.99  # noqa E221
        enforce_margin    =  False                # First and last data point may be excluded in tp  # noqa E221
        use_HCM           =  False                # Use 4 point method, not HCM  # noqa E221
        use_ASTM          =  False                # Use 4 point method, not ASTM  # noqa E221
        residual_method   =  ResidualMethod.NONE  # No processing on residue  # noqa E221  # noqa E221
        spread_damage     =  SDMethod.NONE        # No damage spreading  # noqa E221

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        self.assertEqual(res["rfm"].sum(), 1)
        self.assertEqual(res["rfm"][2 - 1, 3 - 1], 1)
        self.assertTrue((res["res"].flatten() == [4, 1]).all())

    def test_small_sample(self):
        class_count       =  6  # noqa E221
        x                 =  np.array([2, 5, 3, 6, 2, 4, 1, 6, 1,  # noqa E221
                                       4, 1, 5, 3, 6, 3, 6, 1, 5, 2])  # noqa E221
        x_max             =  x.max()  # noqa E221
        x_min             =  x.min()  # noqa E221
        class_width, \
         class_offset     =  self.class_param(x, class_count)  # noqa E221
        hysteresis        =  class_width  # noqa E221
        enforce_margin    =  False                # First and last data point may be excluded in tp  # noqa E221
        use_HCM           =  False                # Use 4 point method, not HCM  # noqa E221
        use_ASTM          =  False                # Use 4 point method, not ASTM  # noqa E221
        residual_method   =  ResidualMethod.NONE  # No processing on residue  # noqa E221
        spread_damage     =  SDMethod.NONE        # No damage spreading  # noqa E221

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        self.assertEqual(res["rfm"].sum(), 7)
        self.assertEqual(res["rfm"][5 - 1, 3 - 1], 2)
        self.assertEqual(res["rfm"][6 - 1, 3 - 1], 1)
        self.assertEqual(res["rfm"][1 - 1, 4 - 1], 1)
        self.assertEqual(res["rfm"][2 - 1, 4 - 1], 1)
        self.assertEqual(res["rfm"][1 - 1, 6 - 1], 2)

        self.assertTrue((res["res"].flatten() == [2, 6, 1, 5, 2]).all())

    def test_long_series(self):
        try:
            import pandas as pd
        except ImportError as err:
            print("This test requires module 'pandas'!")
            raise err

        class_count       =  100  # noqa E221
        class_offset      = -2025  # noqa E221
        class_width       =  50  # noqa E221
        x                 =  pd.read_csv(os.path.join(self.get_script_path(), "long_series.csv"), header=None)  # noqa E221
        x                 =  x.to_numpy().squeeze()  # noqa E221
        hysteresis        =  class_width  # noqa E221
        enforce_margin    =  True                        # First and last data point may be excluded in tp  # noqa E221
        use_HCM           =  False                       # Use 4 point method, not HCM  # noqa E221
        use_ASTM          =  False                       # Use 4 point method, not ASTM  # noqa E221
        residual_method   =  ResidualMethod.NONE         # No processing on residue  # noqa E221
        spread_damage     =  SDMethod.RAMP_AMPLITUDE_23  # No damage spreading  # noqa E221

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        # With residuum:    pd == 9.8934e-06 (repeated)
        # Without residuum: pd == 1.1486e-07
        self.assertTrue(np.absolute(res["tp"][:, 2].sum() / res["damage"] - 1) < 1e-10)  # noqa E501

        spread_damage = SDMethod.TRANSIENT_23c

        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        # assert(abs(sum(dh) / pd - 1) < 1e-10)

        self.assertEqual("%.4e" % res["damage"], "1.1486e-07")
        self.assertEqual(res["rfm"].sum(), 640)
        self.assertEqual(len(res["res"]), 10)
        test = np.absolute(res["res"].flatten() - [
                           0, 142, -609, 2950, -2000,
                           2159, 1894, 2101, 1991, 2061])
        self.assertTrue(test.sum() < 1e-3)


def run():
    unittest.main()
