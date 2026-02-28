"""Unit tests for the rainflow cycle counting implementation.

This module contains test cases for verifying the correctness of the rainflow
cycle counting algorithm, including tests for empty series, single cycles,
small samples, and long data series.
"""

from __future__ import annotations

# ruff: noqa: S101
import logging
import unittest
from pathlib import Path

import numpy as np

from .. import ResidualMethod, SDMethod, rfc

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestRainflowCounting(unittest.TestCase):
    """Unit tests for the rainflow cycle counting implementation.

    This class contains test cases for verifying the correctness of the rainflow
    cycle counting algorithm, including tests for empty series, single cycles,
    small samples, and long data series.
    """

    def get_script_path(self) -> str:
        """Get the directory path of the current script.

        Returns
        -------
        str
            The directory path where the current script is located.

        """
        return Path(__file__).parent

    def class_param(self, data: np.ndarray, class_count: int) -> tuple[float, float]:
        """Calculate class parameters for rainflow cycle counting.

        This method computes the class width and class offset based on the input data
        and the specified number of classes. If the input data is empty, default values
        are used.

        Parameters
        ----------
        data : np.ndarray
            The input array containing the data.
        class_count : int
            The number of classes for the rainflow cycle counting. Must be greater than 1.

        Returns
        -------
        tuple
            A tuple containing the class width and class offset.

        Raises
        ------
        AssertionError
            If `class_count` is not greater than 1.

        """
        if class_count <= 1:
            msg = "class_count must be greater than 1"
            raise ValueError(msg)

        if len(data) == 0:
            class_width = 1  # Default class width when data is empty
            class_offset = 0  # Default class offset when data is empty
        else:
            class_width = np.ptp(data) / (class_count - 1)
            class_width = np.ceil(class_width * 100) / 100  # Round up to the nearest 0.01
            class_offset = np.floor((data.min() - class_width / 2) * 1000) / 1000  # Round down to the nearest 0.001

        return class_width, class_offset

    def test_empty_series(self) -> None:
        """Test the rainflow cycle counting method with an empty data series.

        This test verifies that the rainflow cycle counting method handles an empty
        data series correctly, producing an empty rainflow matrix and no residuals.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If the resulting rainflow matrix sum is not 0 or if residuals are present.

        """
        class_count = 100  # Number of classes
        x = np.array([])  # Empty data array
        x_max = 1  # Placeholder for maximum value in data
        x_min = -1  # Placeholder for minimum value in data

        # Calculate class parameters
        class_width, class_offset = self.class_param(x, class_count)

        hysteresis = class_width  # Hysteresis width
        enforce_margin = False  # Exclude first and last data point in turning points
        use_HCM = False  # Do not use HCM method
        use_ASTM = False  # Do not use ASTM method
        residual_method = ResidualMethod.NONE  # No processing on residue
        spread_damage = SDMethod.NONE  # No damage spreading

        # Perform rainflow counting
        res = rfc(
            x,
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        # Assert that the rainflow matrix sum is 0
        assert res["rfm"].sum() == 0
        # Assert that there are no residuals
        assert len(res["res"]) == 0

    def test_single_cycle_up(self) -> None:
        """Test the rainflow cycle counting method with a single upward cycle.

        This test verifies that the rainflow cycle counting method correctly identifies
        a single cycle in a small dataset. The test checks that the rainflow matrix
        and residuals are accurately computed.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If the resulting rainflow matrix or residuals do not match the expected values.

        """
        class_count = 4  # Number of classes
        x = np.array([1, 3, 2, 4])  # Data array representing a single cycle
        x_max = 4  # Maximum value in data
        x_min = 1  # Minimum value in data

        # Calculate class parameters
        class_width, class_offset = self.class_param(x, class_count)

        hysteresis = class_width * 0.99  # Hysteresis width slightly less than class width
        enforce_margin = False  # Exclude first and last data point in turning points
        use_HCM = False  # Do not use HCM method
        use_ASTM = False  # Do not use ASTM method
        residual_method = ResidualMethod.NONE  # No processing on residue
        spread_damage = SDMethod.NONE  # No damage spreading

        # Perform rainflow counting
        res = rfc(
            x,
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        # Assert that the rainflow matrix sum is 1
        assert res["rfm"].sum() == 1
        # Assert that the specific entry in the rainflow matrix is 1
        assert res["rfm"][3 - 1, 2 - 1] == 1
        # Assert that the residuals match the expected values
        assert (res["res"].flatten() == [1, 4]).all()

    def test_one_cycle_down(self) -> None:
        """Test the rainflow cycle counting method with a single downward cycle.

        This test verifies that the rainflow cycle counting method correctly identifies
        a single cycle in a small dataset with a downward trend. The test checks that the
        rainflow matrix and residuals are accurately computed.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If the resulting rainflow matrix or residuals do not match the expected values.

        """
        class_count = 4  # Number of classes
        x = np.array([4, 2, 3, 1])  # Data array representing a single cycle with a downward trend
        x_max = 4  # Maximum value in data
        x_min = 1  # Minimum value in data

        # Calculate class parameters
        class_width, class_offset = self.class_param(x, class_count)

        hysteresis = class_width * 0.99  # Hysteresis width slightly less than class width
        enforce_margin = False  # Exclude first and last data point in turning points
        use_HCM = False  # Do not use HCM method
        use_ASTM = False  # Do not use ASTM method
        residual_method = ResidualMethod.NONE  # No processing on residue
        spread_damage = SDMethod.NONE  # No damage spreading

        # Perform rainflow counting
        res = rfc(
            x,
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        # Assert that the rainflow matrix sum is 1
        assert res["rfm"].sum() == 1
        # Assert that the specific entry in the rainflow matrix is 1
        assert res["rfm"][2 - 1, 3 - 1] == 1
        # Assert that the residuals match the expected values
        assert (res["res"].flatten() == [4, 1]).all()

    def test_small_sample(self) -> None:
        """Test the rainflow cycle counting method with a small sample dataset.

        This test verifies that the rainflow cycle counting method correctly identifies
        cycles in a small dataset. The test checks that the rainflow matrix and residuals
        are accurately computed.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If the resulting rainflow matrix or residuals do not match the expected values.

        """
        class_count = 6  # Number of classes
        x = np.array([2, 5, 3, 6, 2, 4, 1, 6, 1,
                      4, 1, 5, 3, 6, 3, 6, 1, 5, 2])  # Small sample data array
        x_max = x.max()  # Maximum value in data
        x_min = x.min()  # Minimum value in data

        # Calculate class parameters
        class_width, class_offset = self.class_param(x, class_count)

        hysteresis = class_width  # Hysteresis width
        enforce_margin = False  # Exclude first and last data point in turning points
        use_HCM = False  # Do not use HCM method
        use_ASTM = False  # Do not use ASTM method
        residual_method = ResidualMethod.NONE  # No processing on residue
        spread_damage = SDMethod.NONE  # No damage spreading

        # Perform rainflow counting
        res = rfc(
            x,
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        # Assert that the rainflow matrix sum is 7
        assert res["rfm"].sum() == 7
        # Assert that specific entries in the rainflow matrix are correct
        assert res["rfm"][5 - 1, 3 - 1] == 2
        assert res["rfm"][6 - 1, 3 - 1] == 1
        assert res["rfm"][1 - 1, 4 - 1] == 1
        assert res["rfm"][2 - 1, 4 - 1] == 1
        assert res["rfm"][1 - 1, 6 - 1] == 2

        # Assert that the residuals match the expected values
        assert (res["res"].flatten() == [2, 6, 1, 5, 2]).all()

    def test_long_series(self) -> None:
        """Test the rainflow cycle counting method with a long data series.

        This test verifies that the rainflow cycle counting method correctly processes
        a long dataset read from a CSV file. The test checks the computed damage, rainflow
        matrix, and residuals for accuracy.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        ImportError
            If the required module 'pandas' is not installed.
        AssertionError
            If the resulting damage, rainflow matrix, or residuals do not match the expected values.

        """
        try:
            import pandas as pd
        except ImportError:
            logger.exception("This test requires module 'pandas'!")
            raise

        class_count = 100  # Number of classes
        class_offset = -2025  # Class offset
        class_width = 50  # Class width

        # Read the data from CSV
        x = pd.read_csv(Path(self.get_script_path()) / "long_series.csv", header=None)
        x = x.to_numpy().squeeze()  # Convert to numpy array and squeeze

        hysteresis = class_width  # Hysteresis width
        enforce_margin = True  # Exclude first and last data point in turning points
        use_HCM = False  # Do not use HCM method
        use_ASTM = False  # Do not use ASTM method
        residual_method = ResidualMethod.NONE  # No processing on residue
        spread_damage = SDMethod.RAMP_AMPLITUDE_23  # Spread damage method

        # Perform rainflow counting
        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        # With residuum:    pd == 9.8934e-06 (repeated)
        # Without residuum: pd == 1.1486e-07
        assert np.absolute(res["tp"][:, 2].sum() / res["damage"] - 1) < 1e-10

        # Change spread damage method
        spread_damage = SDMethod.TRANSIENT_23c

        # Perform rainflow counting again with the new spread damage method
        res = rfc(
            x, class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=hysteresis,
            residual_method=residual_method,
            enforce_margin=enforce_margin,
            use_HCM=use_HCM,
            use_ASTM=use_ASTM,
            spread_damage=spread_damage,
        )

        assert "{:.4e}".format(res["damage"]) == "1.1486e-07"
        assert res["rfm"].sum() == 640
        assert len(res["res"]) == 10

        test = np.absolute(res["res"].flatten() - [
            0, 142, -609, 2950, -2000,
            2159, 1894, 2101, 1991, 2061,
        ])
        assert test.sum() < 1e-3


def run() -> None:
    """Run all unit tests in this module."""
    unittest.main()
