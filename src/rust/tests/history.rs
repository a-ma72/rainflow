mod common;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use common::{assert_matrix, assert_near, assert_vec};
use rfcnt::{rfc, Options, ResidualMethod, Rfc, SdMethod, Woehler};

#[test]
fn spread_methods_sum_to_damage_and_differ() {
    let mut data = Vec::new();
    for value in 0..=100 {
        data.push(value as f64);
    }
    for value in (0..100).rev() {
        data.push(value as f64);
    }
    for value in 0..=40 {
        data.push(value as f64);
    }
    for value in (0..40).rev() {
        data.push(value as f64);
    }
    for value in 0..=80 {
        data.push(value as f64);
    }

    let mut opt = Options::default();
    opt.class_count = 120;
    opt.class_offset = -0.5;
    opt.hysteresis = Some(1.0);
    opt.enforce_margin = false;
    // Residue handling stays off. Repeated counting indexes past the input buffer.
    opt.residual_method = ResidualMethod::None;
    opt.wl = Woehler::Modified { sx: 50.0, nx: 1e3, k: 3.0, k2: 3.0 };
    opt.spread_damage = SdMethod::Transient23c;
    let transient = rfc(&data, 1.0, &opt).unwrap();

    opt.spread_damage = SdMethod::RampAmplitude23;
    let ramp = rfc(&data, 1.0, &opt).unwrap();

    for res in [&transient, &ramp] {
        assert_eq!(res.dh.len(), data.len());
        let sum: f64 = res.dh.iter().sum();
        assert_near(sum, res.damage);
        assert!(res.damage > 0.0, "damage {}", res.damage);
    }
    let max_gap = transient
        .dh
        .iter()
        .zip(&ramp.dh)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        max_gap > 0.01 * transient.damage,
        "gap {max_gap} transient {} ramp {}",
        transient.damage,
        ramp.damage
    );
}

fn load_long_series() -> Vec<f64> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../python/tests/long_series.csv");
    let file = File::open(&path).unwrap_or_else(|err| panic!("open {}: {err}", path.display()));
    BufReader::new(file)
        .lines()
        .map(|line| {
            let line = line.expect("read long series");
            line.trim().parse::<f64>().unwrap_or_else(|err| panic!("parse {line:?}: {err}"))
        })
        .collect()
}

fn long_options(spread: SdMethod, residual: ResidualMethod) -> Options {
    let mut opt = Options::default();
    opt.class_count = 100;
    opt.class_offset = -2025.0;
    opt.hysteresis = Some(50.0);
    opt.enforce_margin = true;
    opt.residual_method = residual;
    opt.spread_damage = spread;
    opt
}

#[test]
fn long_series_matches_python_reference() {
    let data = load_long_series();
    let width = 50.0;

    let ramp = rfc(&data, width, &long_options(SdMethod::RampAmplitude23, ResidualMethod::None)).unwrap();
    let tp_damage: f64 = ramp.tp.iter().map(|point| point.damage).sum();
    assert!((tp_damage / ramp.damage - 1.0).abs() < 1e-10);

    let res = rfc(&data, width, &long_options(SdMethod::Transient23c, ResidualMethod::None)).unwrap();
    assert!(
        (res.damage - 1.1486e-7).abs() <= 5e-12,
        "damage {:.10e}",
        res.damage
    );
    assert_near(res.rfm.sum(), 640.0);
    assert_eq!(res.res.len(), 10);
    let expect = [0.0, 142.0, -609.0, 2950.0, -2000.0, 2159.0, 1894.0, 2101.0, 1991.0, 2061.0];
    let delta: f64 = res.res.iter().zip(expect).map(|(got, exp)| (got - exp).abs()).sum();
    assert!(delta < 1e-3, "residue delta {delta}");
}

#[test]
fn long_series_chunked_feed_matches_one_shot() {
    let data = load_long_series();
    let width = 50.0;
    let opt = long_options(SdMethod::None, ResidualMethod::None);
    let one = rfc(&data, width, &opt).unwrap();

    let mut rf = Rfc::new(width, &opt).unwrap();
    for chunk in data.chunks(1000) {
        rf.feed(chunk).unwrap();
    }
    rf.finalize(ResidualMethod::None).unwrap();

    assert_near(rf.damage().unwrap(), one.damage);
    assert_matrix(&rf.rfm().unwrap(), &one.rfm);
    assert_vec(&rf.residue().unwrap(), &one.res);
    assert_vec(&rf.res_raw().unwrap(), &one.res_raw);
}
