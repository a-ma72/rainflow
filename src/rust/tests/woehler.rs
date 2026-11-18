mod common;

use common::{assert_near, curve_series};
use rfcnt::{damage_from_rp, rfc, Options, ResidualMethod, RpDamageCalcMethod, Woehler};

fn base_options() -> Options {
    let mut opt = Options::for_counter();
    opt.class_count = 10;
    opt.class_offset = -0.5;
    opt.hysteresis = Some(1.0);
    opt.residual_method = ResidualMethod::Repeated;
    opt
}

#[test]
fn steeper_slope_changes_damage_and_is_reported() {
    let mut opt = base_options();
    let baseline = rfc(&curve_series(), 1.0, &opt).unwrap();
    opt.wl = Woehler::Modified { sx: 1e3, nx: 1e7, k: 8.0, k2: 8.0 };
    let steep = rfc(&curve_series(), 1.0, &opt).unwrap();
    assert!((baseline.damage - steep.damage).abs() > 0.0);
    assert_near(steep.wl_miner_consistent.k.abs(), 8.0);
}

#[test]
fn omission_drops_the_small_cycle() {
    let mut opt = base_options();
    opt.wl = Woehler::Full {
        sx: 1e3,
        nx: 1e7,
        k: 5.0,
        k2: 5.0,
        sd: 0.0,
        nd: f64::MAX,
        omission: 0.0,
    };
    let kept = rfc(&curve_series(), 1.0, &opt).unwrap();
    opt.wl = Woehler::Full {
        sx: 1e3,
        nx: 1e7,
        k: 5.0,
        k2: 5.0,
        sd: 0.0,
        nd: f64::MAX,
        omission: 1.5,
    };
    let omitted = rfc(&curve_series(), 1.0, &opt).unwrap();
    assert!(omitted.damage < kept.damage);
    assert!(kept.damage > 0.0);
}

#[test]
fn damage_from_rp_sorts_and_distinguishes_miner_rules() {
    let wl = Woehler::Full {
        sx: 100.0,
        nx: 1e6,
        k: 5.0,
        k2: 9.0,
        sd: 20.0,
        nd: 1e7,
        omission: 0.0,
    };
    let sa = [5.0, 15.0, 40.0, 80.0, 150.0];
    let counts = [1.0e5, 1.0e4, 1000.0, 100.0, 10.0];
    let sorted = damage_from_rp(&sa, &counts, Some(&wl), RpDamageCalcMethod::MinerElementar).unwrap();
    let shuffled = damage_from_rp(
        &[150.0, 5.0, 80.0, 15.0, 40.0],
        &[10.0, 1.0e5, 100.0, 1.0e4, 1000.0],
        Some(&wl),
        RpDamageCalcMethod::MinerElementar,
    )
    .unwrap();
    assert_near(sorted, shuffled);

    let methods = [
        RpDamageCalcMethod::Default,
        RpDamageCalcMethod::MinerElementar,
        RpDamageCalcMethod::MinerModified,
        RpDamageCalcMethod::MinerConsistent,
    ];
    let mut values = Vec::new();
    for method in methods {
        let damage = damage_from_rp(&sa, &counts, Some(&wl), method).unwrap();
        assert!(damage.is_finite() && damage > 0.0, "method {} -> {damage}", method as i32);
        values.push(damage);
    }
    assert!(values.iter().any(|value| !near_loose(*value, values[0])));
}

fn near_loose(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-6 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn damage_from_rp_rejects_histogram_above_class_limit() {
    let sa = vec![1.0; 1025];
    let counts = vec![1.0; 1025];
    assert!(damage_from_rp(&sa, &counts, None, RpDamageCalcMethod::Default).is_err());
}
