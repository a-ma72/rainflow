import os
import unittest
from rfcnt import rfcnt
import numpy as np


class TestRainflowCounting(unittest.TestCase):

    def get_script_path(self):
        return os.path.dirname(__file__)

    def class_param(self, data, class_count):
        assert(class_count > 1)

        if len(data) == 0:
            class_width  = 1
            class_offset = 0
        else:
            class_width  = (data.max() - data.min()) / (class_count - 1)
            class_width  = np.ceil(class_width * 100) / 100
            class_offset = np.floor((data.min() - class_width / 2) * 1000) / 1000

        return class_width, class_offset

    def test_empty_series(self):
        class_count       =  100
        x                 =  np.array([])
        x_max             =  1
        x_min             = -1
        class_width, \
         class_offset     =  self.class_param(x, class_count)
        hysteresis        =  class_width
        enforce_margin    =  0  # First and last data point may be excluded in tp
        use_HCM           =  0  # Use 4 point method, not HCM
        use_ASTM          =  0  # Use 4 point method, not ASTM
        residual_method   =  0  # No processing on residue
        spread_damage     =  0  # No damage spreading

        res = rfcnt.rfc(
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
        class_count       =  4
        x                 =  np.array([1, 3, 2, 4])
        x_max             =  4
        x_min             =  1
        class_width, \
         class_offset     =  self.class_param(x, class_count)
        hysteresis        =  class_width * 0.99
        enforce_margin    =  0  # First and last data point may be excluded in tp
        use_HCM           =  0  # Use 4 point method, not HCM
        use_ASTM          =  0  # Use 4 point method, not ASTM
        residual_method   =  0  # No processing on residue
        spread_damage     =  0  # No damage spreading

        res = rfcnt.rfc(
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
        class_count       =  4
        x                 =  np.array([4, 2, 3, 1])
        x_max             =  4
        x_min             =  1
        class_width, \
         class_offset     =  self.class_param(x, class_count)
        hysteresis        =  class_width * 0.99
        enforce_margin    =  0  # First and last data point may be excluded in tp
        use_HCM           =  0  # Use 4 point method, not HCM
        use_ASTM          =  0  # Use 4 point method, not ASTM
        residual_method   =  0  # No processing on residue
        spread_damage     =  0  # No damage spreading

        res = rfcnt.rfc(
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
        class_count       =  6
        x                 =  np.array([2, 5, 3, 6, 2, 4, 1, 6, 1,
                                       4, 1, 5, 3, 6, 3, 6, 1, 5, 2])
        x_max             =  x.max()
        x_min             =  x.min()
        class_width, \
         class_offset     =  self.class_param(x, class_count)
        hysteresis        =  class_width
        enforce_margin    =  0  # First and last data point may be excluded in tp
        use_HCM           =  0  # Use 4 point method, not HCM
        use_ASTM          =  0  # Use 4 point method, not ASTM
        residual_method   =  0  # No processing on residue
        spread_damage     =  0  # No damage spreading

        res = rfcnt.rfc(
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

        class_count       =  100
        class_offset      = -2025
        class_width       =  50
        x                 =  pd.read_csv(os.path.join(self.get_script_path(), "long_series.csv"), header=None)
        x                 =  x.to_numpy().squeeze()
        hysteresis        =  class_width
        use_HCM           =  0
        enforce_margin    =  1  # Enforce first and last data point included in tp
        use_HCM           =  0  # Use 4 point method, not HCM
        use_ASTM          =  0  # Use 4 point method, not ASTM
        residual_method   =  0  # 0=RFC_RES_NONE, 7=RFC_RES_REPEATED
        spread_damage     =  1  # 0=RFC_SD_HALF_23, 1=RFC_SD_RAMP_AMPLITUDE_23

        res = rfcnt.rfc(
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
        self.assertTrue(np.absolute(res["tp"][:, 2].sum() / res["damage"] - 1) < 1e-10)

        spread_damage = 8;  # 7=RFC_SD_TRANSIENT_23, 8=RFC_SD_TRANSIENT_23c

        res = rfcnt.rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage)

        #assert(abs(sum(dh) / pd - 1) < 1e-10)

        self.assertEqual("%.4e" % res["damage"], "1.1486e-07")
        self.assertEqual(res["rfm"].sum(), 640)
        self.assertEqual(len(res["res"]), 10)
        test = np.absolute(res["res"].flatten() - [
                           0, 142, -609, 2950, -2000,
                           2159, 1894, 2101, 1991, 2061])
        self.assertTrue(test.sum() < 1e-3)


def run():
    unittest.main()
