mod common;

use common::{assert_near, curve_series};
use rfcnt::{rfc, Options, ResidualMethod, RfcResult, Woehler};

fn cycle_damage(from: f64, to: f64, sx: f64, nx: f64, k: f64) -> f64 {
    let amplitude = (from - to).abs() / 2.0;
    (amplitude / sx).powf(k.abs()) / nx
}

fn count(method: ResidualMethod) -> RfcResult {
    let mut opt = Options::for_counter();
    opt.class_count = 10;
    opt.class_offset = -0.5;
    opt.hysteresis = Some(1.0);
    opt.residual_method = method;
    opt.wl = Woehler::Modified { sx: 1e3, nx: 1e7, k: 5.0, k2: 5.0 };
    rfc(&curve_series(), 1.0, &opt).unwrap()
}

#[test]
fn residue_methods_on_inner_and_outer_cycle() {
    let sx = 1e3;
    let nx = 1e7;
    let k = 5.0;
    let inner = cycle_damage(5.0, 3.0, sx, nx, k);
    let outer = cycle_damage(7.0, 2.0, sx, nx, k);

    let repeated = count(ResidualMethod::Repeated);
    let half = count(ResidualMethod::Halfcycles);
    let full = count(ResidualMethod::Fullcycles);
    let ignore = count(ResidualMethod::Ignore);

    assert!(repeated.res.is_empty());
    assert_near(repeated.damage, inner + outer);

    // The inner cycle is closed in full; the outer residue is counted as a half cycle.
    assert_near(half.damage, inner + outer / 2.0);
    assert!(half.damage < repeated.damage);

    assert_near(full.damage, inner + outer);
    assert!(full.damage > half.damage);

    assert!(ignore.damage < repeated.damage);
    assert!(!ignore.res.is_empty());
}

#[test]
fn din45667_range_pair_counts() {
    let data = [4.9, 6.0, 4.0, 7.0, 3.0, 9.0, 5.0, 8.0, 6.9];
    let mut opt = Options::for_counter();
    opt.class_count = 10;
    opt.class_offset = -0.5;
    opt.hysteresis = Some(1.0);
    opt.residual_method = ResidualMethod::Din45667;
    opt.wl = Woehler::Modified { sx: 1e3, nx: 1e7, k: 5.0, k2: 5.0 };

    let res = rfc(&data, 1.0, &opt).unwrap();
    let counts: Vec<f64> = res.rp.iter().map(|row| row[1]).collect();
    assert_eq!(counts.len(), 10);
    let expect = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for (got, expect) in counts.iter().zip(expect) {
        assert_near(*got, expect);
    }
}
