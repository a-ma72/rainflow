"""Unit tests for the stateful RFC class."""

from __future__ import annotations

# ruff: noqa: S101
import unittest

import numpy as np

import rfcnt

RFC = rfcnt.RFC
rfc = rfcnt.rfc

try:
    ResidualMethod = rfcnt.ResidualMethod
    SDMethod = rfcnt.SDMethod
except AttributeError:
    class ResidualMethod:
        NONE = 0
        REPEATED = 7

    class SDMethod:
        NONE = -1


class TestRFCClass(unittest.TestCase):
    """Tests for incremental feeding via rfcnt.RFC."""

    @staticmethod
    def _class_param(data: np.ndarray, class_count: int) -> tuple[float, float]:
        class_width = np.ptp(data) / (class_count - 1)
        class_width = np.ceil(class_width * 100) / 100
        class_offset = np.floor((data.min() - class_width / 2) * 1000) / 1000
        return float(class_width), float(class_offset)

    @staticmethod
    def _assert_wl_miner_consistent(rf: object, one_shot: dict) -> None:
        got = dict(rf.wl_miner_consistent)  # type: ignore[attr-defined]
        expected = dict(one_shot["wl_miner_consistent"])
        if set(got) != set(expected):
            raise AssertionError(f"wl keys {set(got)} != {set(expected)}")
        for key, value in expected.items():
            np.testing.assert_allclose(got[key], value, equal_nan=True)

    def test_multi_feed_matches_one_shot_rfc(self) -> None:
        """Chunked feed + finalize matches one-shot rfc() with DH off."""
        class_count = 4
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0])
        class_width, class_offset = self._class_param(x, class_count)
        kwargs = dict(
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
        )

        one_shot = rfc(x, **kwargs)

        rf = RFC(
            class_width,
            class_count=class_count,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
        )
        split = len(x) // 2
        rf.feed(x[:split])
        rf.feed(x[split:])
        rf.finalize(ResidualMethod.REPEATED)

        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.residue).flatten(), np.asarray(one_shot["res"]).flatten())
        np.testing.assert_allclose(np.asarray(rf.rp), one_shot["rp"])
        np.testing.assert_allclose(np.asarray(rf.lc), one_shot["lc"])
        np.testing.assert_allclose(np.asarray(rf.rfm), one_shot["rfm"])
        np.testing.assert_allclose(np.asarray(rf.tp), one_shot["tp"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), np.asarray(one_shot["res_raw"]).flatten())
        self._assert_wl_miner_consistent(rf, one_shot)

    def test_feed_after_finalize_raises(self) -> None:
        x = np.array([0.0, 1.0, 0.5, 2.0])
        class_width, class_offset = self._class_param(x, 10)
        rf = RFC(
            class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        )
        rf.feed(x)
        rf.finalize()
        with self.assertRaises(RuntimeError):
            rf.feed(np.array([1.0]))

    def test_residue_readable_before_finalize(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0])
        class_width, class_offset = self._class_param(x, 4)
        rf = RFC(
            class_width,
            class_count=4,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        )
        rf.feed(x[:3])
        residue = np.asarray(rf.residue).flatten()
        self.assertGreaterEqual(residue.size, 1)
        d_before = rf.damage
        rf.feed(x[3:])
        self.assertGreaterEqual(rf.damage, d_before)
        rf.finalize(ResidualMethod.REPEATED)

    def test_damage_as_matches_later_finalize(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0])
        class_width, class_offset = self._class_param(x, 4)
        rf = RFC(
            class_width,
            class_count=4,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
        )
        split = len(x) // 2
        rf.feed(x[:split])
        rf.feed(x[split:])
        preview = rf.damage_as(ResidualMethod.REPEATED)
        rp_preview = np.asarray(rf.rp_as(ResidualMethod.REPEATED))
        lc_preview = np.asarray(rf.lc_as(ResidualMethod.REPEATED))
        rfm_preview = np.asarray(rf.rfm_as(ResidualMethod.REPEATED))
        one_shot = rfc(
            x,
            class_count=4,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
        )
        rf.finalize(ResidualMethod.REPEATED)
        self.assertAlmostEqual(rf.damage, preview)
        np.testing.assert_allclose(np.asarray(rf.rp), rp_preview)
        np.testing.assert_allclose(np.asarray(rf.lc), lc_preview)
        np.testing.assert_allclose(np.asarray(rf.rfm), rfm_preview)
        self.assertAlmostEqual(preview, one_shot["damage"])
        np.testing.assert_allclose(rp_preview, one_shot["rp"])
        np.testing.assert_allclose(lc_preview, one_shot["lc"])
        np.testing.assert_allclose(rfm_preview, one_shot["rfm"])

    def test_damage_as_does_not_mutate(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0])
        class_width, class_offset = self._class_param(x, 4)
        rf = RFC(
            class_width,
            class_count=4,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        )
        rf.feed(x[:4])
        state_before = rf.state
        residue_before = np.asarray(rf.residue).flatten().copy()
        live_before = rf.damage
        rp_before = np.asarray(rf.rp).copy()
        rfm_before = np.asarray(rf.rfm).copy()
        _ = rf.damage_as(ResidualMethod.REPEATED)
        _ = rf.rp_as(ResidualMethod.REPEATED)
        _ = rf.rfm_as(ResidualMethod.REPEATED)
        self.assertEqual(rf.state, state_before)
        np.testing.assert_allclose(np.asarray(rf.residue).flatten(), residue_before)
        self.assertAlmostEqual(rf.damage, live_before)
        np.testing.assert_allclose(np.asarray(rf.rp), rp_before)
        np.testing.assert_allclose(np.asarray(rf.rfm), rfm_before)
        rf.feed(x[4:])
        rf.finalize(ResidualMethod.REPEATED)

    def test_damage_as_after_finalize_raises(self) -> None:
        x = np.array([0.0, 1.0, 0.5, 2.0])
        class_width, class_offset = self._class_param(x, 10)
        rf = RFC(
            class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        )
        rf.feed(x)
        rf.finalize()
        with self.assertRaises(RuntimeError):
            rf.damage_as()
        with self.assertRaises(RuntimeError):
            rf.rp_as()

        rf2 = RFC(
            class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        )
        rf2.feed(x)
        rf2.close()
        with self.assertRaises(RuntimeError):
            rf2.damage_as()
        with self.assertRaises(RuntimeError):
            rf2.rp_as()

    def test_hcm_and_astm_mutually_exclusive(self) -> None:
        with self.assertRaises(ValueError):
            RFC(1.0, class_count=4, use_HCM=True, use_ASTM=True)

    def _assert_rfc_matches_one_shot(self, x: np.ndarray, extra: dict) -> None:
        class_count = extra.get("class_count", 4)
        class_width, class_offset = self._class_param(x, class_count)
        kwargs = dict(
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            **{k: v for k, v in extra.items() if k != "class_count"},
        )
        one_shot = rfc(x, **kwargs)
        rfc_only = {"residual_method", "spread_damage"}
        rf = RFC(**{k: v for k, v in kwargs.items() if k not in rfc_only})
        split = max(len(x) // 2, 1)
        rf.feed(x[:split])
        rf.feed(x[split:])
        preview = rf.damage_as(ResidualMethod.REPEATED)
        res_raw_before = np.asarray(rf.res_raw).flatten().copy()
        rf.finalize(ResidualMethod.REPEATED)
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), res_raw_before)
        self.assertAlmostEqual(preview, one_shot["damage"])
        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.rp), one_shot["rp"])
        np.testing.assert_allclose(np.asarray(rf.lc), one_shot["lc"])
        np.testing.assert_allclose(np.asarray(rf.rfm), one_shot["rfm"])
        np.testing.assert_allclose(np.asarray(rf.tp), one_shot["tp"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), np.asarray(one_shot["res_raw"]).flatten())
        self._assert_wl_miner_consistent(rf, one_shot)

    def test_hcm_matches_one_shot_rfc(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0])
        self._assert_rfc_matches_one_shot(x, {"use_HCM": True})

    def test_astm_matches_one_shot_rfc(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0])
        self._assert_rfc_matches_one_shot(x, {"use_ASTM": True})

    def test_res_raw_stable_across_finalize(self) -> None:
        x = np.array([1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0])
        class_width, class_offset = self._class_param(x, 4)
        rf = RFC(
            class_width,
            class_count=4,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
        )
        rf.feed(x)
        before = np.asarray(rf.res_raw).flatten().copy()
        self.assertTrue(np.asarray(rf.res_raw).flags.writeable)
        rf.finalize(ResidualMethod.REPEATED)
        snap = np.asarray(rf.res_raw)
        np.testing.assert_allclose(snap.flatten(), before)
        self.assertFalse(snap.flags.writeable)
        if snap.size:
            with self.assertRaises(ValueError):
                snap.reshape(-1)[0] = 0.0
            snap.flags.writeable = True
            snap.reshape(-1)[0] = 0.0
            np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), before)

    def test_hcm_autoresize_matches_one_shot_and_damage_as(self) -> None:
        x = np.array([1.0, 3.0, 0.5, 10.0, 2.0, 8.0, 1.0])
        class_count = 4
        class_width, class_offset = self._class_param(x[:3], class_count)
        kwargs = dict(
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
            use_HCM=True,
        )
        one_shot = rfc(x, **kwargs)
        rf = RFC(
            class_width,
            class_count=class_count,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
            use_HCM=True,
        )
        rf.feed(x[:3])
        rf.feed(x[3:])
        preview = rf.damage_as(ResidualMethod.REPEATED)
        res_raw_before = np.asarray(rf.res_raw).flatten().copy()
        rf.finalize(ResidualMethod.REPEATED)
        self.assertAlmostEqual(preview, rf.damage)
        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), res_raw_before)
        np.testing.assert_allclose(np.asarray(rf.rp), one_shot["rp"])
        np.testing.assert_allclose(np.asarray(rf.rfm), one_shot["rfm"])
        np.testing.assert_allclose(np.asarray(rf.tp), one_shot["tp"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), np.asarray(one_shot["res_raw"]).flatten())
        self._assert_wl_miner_consistent(rf, one_shot)

    def test_astm_autoresize_matches_one_shot_and_damage_as(self) -> None:
        x = np.array([1.0, 3.0, 0.5, 10.0, 2.0, 8.0, 1.0])
        class_count = 4
        class_width, class_offset = self._class_param(x[:3], class_count)
        kwargs = dict(
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
            use_ASTM=True,
        )
        one_shot = rfc(x, **kwargs)
        rf = RFC(
            class_width,
            class_count=class_count,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
            use_ASTM=True,
        )
        rf.feed(x[:3])
        rf.feed(x[3:])
        preview = rf.damage_as(ResidualMethod.REPEATED)
        res_raw_before = np.asarray(rf.res_raw).flatten().copy()
        rf.finalize(ResidualMethod.REPEATED)
        self.assertAlmostEqual(preview, rf.damage)
        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), res_raw_before)
        np.testing.assert_allclose(np.asarray(rf.rp), one_shot["rp"])
        np.testing.assert_allclose(np.asarray(rf.rfm), one_shot["rfm"])
        np.testing.assert_allclose(np.asarray(rf.tp), one_shot["tp"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), np.asarray(one_shot["res_raw"]).flatten())
        self._assert_wl_miner_consistent(rf, one_shot)

    def test_autoresize_matches_one_shot_and_damage_as(self) -> None:
        x = np.array([1.0, 3.0, 0.5, 10.0, 2.0, 8.0, 1.0])
        class_count = 4
        class_width, class_offset = self._class_param(x[:3], class_count)
        kwargs = dict(
            class_count=class_count,
            class_width=class_width,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            residual_method=ResidualMethod.REPEATED,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
        )
        one_shot = rfc(x, **kwargs)
        rf = RFC(
            class_width,
            class_count=class_count,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            enforce_margin=False,
            spread_damage=SDMethod.NONE,
            auto_resize=True,
        )
        rf.feed(x[:3])
        rf.feed(x[3:])
        preview = rf.damage_as(ResidualMethod.REPEATED)
        rf.finalize(ResidualMethod.REPEATED)
        self.assertAlmostEqual(preview, rf.damage)
        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.rp), one_shot["rp"])
        np.testing.assert_allclose(np.asarray(rf.rfm), one_shot["rfm"])
        np.testing.assert_allclose(np.asarray(rf.tp), one_shot["tp"])
        np.testing.assert_allclose(np.asarray(rf.res_raw).flatten(), np.asarray(one_shot["res_raw"]).flatten())
        self._assert_wl_miner_consistent(rf, one_shot)

    def test_spread_damage_rejected(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            RFC(1.0, class_count=4, spread_damage=SDMethod.TRANSIENT_23c)
        self.assertIn("rfc()", str(ctx.exception))
        with self.assertRaises(TypeError):
            RFC(1.0, class_count=4, dh_capacity=8)

    def test_lc_method_slopes_all(self) -> None:
        x = np.array([0.0, 1.0, 0.0, 2.0])
        class_width, class_offset = self._class_param(x, 10)
        kwargs = dict(
            class_width=class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
            lc_method=rfcnt.LCMethod.SLOPES_ALL,
        )
        one_shot = rfc(x, residual_method=ResidualMethod.REPEATED, enforce_margin=False, **kwargs)
        rf = RFC(enforce_margin=False, **kwargs)
        rf.feed(x)
        rf.finalize(ResidualMethod.REPEATED)
        self.assertAlmostEqual(rf.damage, one_shot["damage"])
        np.testing.assert_allclose(np.asarray(rf.lc), one_shot["lc"])

    def test_lc_method_default_matches_slopes_all(self) -> None:
        x = np.array([0.0, 1.0, 0.0, 2.0])
        class_width, class_offset = self._class_param(x, 10)
        kwargs = dict(
            class_width=class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
            enforce_margin=False,
        )
        default = rfc(x, residual_method=ResidualMethod.REPEATED, **kwargs)
        both = rfc(
            x,
            residual_method=ResidualMethod.REPEATED,
            lc_method=rfcnt.LCMethod.SLOPES_ALL,
            **kwargs,
        )
        up = rfc(
            x,
            residual_method=ResidualMethod.REPEATED,
            lc_method=rfcnt.LCMethod.SLOPES_UP,
            **kwargs,
        )
        np.testing.assert_allclose(default["lc"], both["lc"])
        self.assertFalse(np.allclose(default["lc"], up["lc"]))

    def test_lc_method_rejects_out_of_range(self) -> None:
        """Values outside LCMethod 0..3 are rejected (3 is FVA, not a mask)."""
        x = np.array([0.0, 1.0, 0.0, 2.0])
        with self.assertRaises(ValueError) as ctx:
            RFC(1.0, class_count=4, lc_method=4)
        self.assertIn("0, 1, 2 or 3", str(ctx.exception))
        with self.assertRaises(ValueError):
            rfc(x, class_width=1.0, class_count=4, lc_method=4)

    def test_lc_method_fva_t1(self) -> None:
        x = np.array([-3.0, 2.0, -1.0])
        kwargs = dict(
            class_width=1.0,
            class_count=6,
            class_offset=-3.5,
            hysteresis=0.5,
            spread_damage=SDMethod.NONE,
            lc_method=rfcnt.LCMethod.FVA,
            enforce_margin=True,
        )
        expect = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
        expect_n = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 0.0])
        one_shot = rfc(x, residual_method=ResidualMethod.REPEATED, **kwargs)
        np.testing.assert_allclose(one_shot["lc"][:, 0], expect)
        np.testing.assert_allclose(one_shot["lc"][:, 1], expect_n)

        rf = RFC(**kwargs)
        rf.feed(x)
        np.testing.assert_allclose(np.asarray(rf.lc)[:, 1], expect_n)
        np.testing.assert_allclose(np.asarray(rf.lc_as(ResidualMethod.REPEATED))[:, 1], expect_n)
        rf.finalize(ResidualMethod.REPEATED)
        np.testing.assert_allclose(np.asarray(rf.lc)[:, 1], expect_n)

    def test_lc_method_fva_t2(self) -> None:
        x = np.array([1.0, -2.0, 3.0])
        kwargs = dict(
            class_width=1.0,
            class_count=6,
            class_offset=-2.5,
            hysteresis=0.5,
            spread_damage=SDMethod.NONE,
            lc_method=rfcnt.LCMethod.FVA,
            enforce_margin=True,
        )
        expect_n = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
        one_shot = rfc(x, residual_method=ResidualMethod.REPEATED, **kwargs)
        np.testing.assert_allclose(one_shot["lc"][:, 1], expect_n)

        rf = RFC(**kwargs)
        rf.feed(x)
        np.testing.assert_allclose(np.asarray(rf.lc)[:, 1], expect_n)
        np.testing.assert_allclose(np.asarray(rf.lc_as(ResidualMethod.NONE))[:, 1], expect_n)
        rf.finalize(ResidualMethod.REPEATED)
        np.testing.assert_allclose(np.asarray(rf.lc)[:, 1], expect_n)

    def test_lc_method_din45667_is_fva_alias(self) -> None:
        self.assertEqual(rfcnt.LCMethod.DIN45667, rfcnt.LCMethod.FVA)
        self.assertEqual(int(rfcnt.LCMethod.DIN45667), 3)

    def test_lc_method_fva_sign_dependent_zero_crossing(self) -> None:
        """FVA switches direction at zero; DIN SLOPES_UP does not."""
        x = np.array([-3.0, 2.0, -1.0, 2.5, -2.0])
        kwargs = dict(
            class_width=1.0,
            class_count=7,
            class_offset=-3.5,
            hysteresis=0.5,
            spread_damage=SDMethod.NONE,
            enforce_margin=True,
        )
        up = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.SLOPES_UP, **kwargs)
        dn = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.SLOPES_DOWN, **kwargs)
        fva = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.FVA, **kwargs)
        levels = fva["lc"][:, 0]
        expect = np.where(levels >= 0.0, up["lc"][:, 1], dn["lc"][:, 1])
        np.testing.assert_allclose(fva["lc"][:, 1], expect)
        self.assertFalse(np.allclose(up["lc"][:, 1], fva["lc"][:, 1]))

    def test_lc_method_fva_zero_class_bound(self) -> None:
        """Exact u=0 follows the positive (rising) FVA branch."""
        x = np.array([-2.0, 2.0, -2.0])
        kwargs = dict(
            class_width=1.0,
            class_count=5,
            class_offset=-2.0,
            hysteresis=0.5,
            spread_damage=SDMethod.NONE,
            enforce_margin=True,
        )
        up = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.SLOPES_UP, **kwargs)
        dn = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.SLOPES_DOWN, **kwargs)
        fva = rfc(x, residual_method=ResidualMethod.REPEATED, lc_method=rfcnt.LCMethod.FVA, **kwargs)
        levels = fva["lc"][:, 0]
        zero = np.isclose(levels, 0.0)
        self.assertTrue(np.any(zero))
        np.testing.assert_allclose(fva["lc"][zero, 1], up["lc"][zero, 1])
        expect = np.where(levels >= 0.0, up["lc"][:, 1], dn["lc"][:, 1])
        np.testing.assert_allclose(fva["lc"][:, 1], expect)

    def test_data_must_be_1d(self) -> None:
        x2 = np.array([[0.0, 1.0], [0.0, 2.0]])
        with self.assertRaises(ValueError) as ctx:
            rfc(x2, class_width=1.0, class_count=4)
        self.assertIn("1-D", str(ctx.exception))
        rf = RFC(1.0, class_count=4)
        with self.assertRaises(ValueError):
            rf.feed(x2)

    def test_context_manager_closes(self) -> None:
        x = np.array([0.0, 1.0, 0.0, 2.0])
        class_width, class_offset = self._class_param(x, 10)
        with RFC(
            class_width,
            class_count=10,
            class_offset=class_offset,
            hysteresis=class_width * 0.99,
            spread_damage=SDMethod.NONE,
        ) as rf:
            rf.feed(x)
            self.assertGreaterEqual(rf.state, 1)
        with self.assertRaises(RuntimeError):
            _ = rf.damage

    def test_at_init_transform(self) -> None:
        rf = RFC(1.0, class_count=10, class_offset=0.0, hysteresis=1.0,
                 spread_damage=SDMethod.NONE)
        rf.at_init(0.3, R_rig=-1.0, R_pinned=True)
        got = rf.at_transform(np.array([3.0, 2.0]), np.array([1.0, 2.0]))
        np.testing.assert_allclose(got, [3.3, 2.6], atol=1e-10)

        rf.feed(np.array([0.0, 1.0, 0.0]))
        with self.assertRaises(RuntimeError):
            rf.at_init(0.3)


if __name__ == "__main__":
    unittest.main()
