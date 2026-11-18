#![allow(dead_code)]

use rfcnt::{Matrix, TurningPoint};

/// Class width and offset, matching the Python `class_param` helper.
pub fn class_param(data: &[f64], class_count: u32) -> (f64, f64) {
    if data.is_empty() {
        return (1.0, 0.0);
    }
    let min = data.iter().copied().fold(f64::INFINITY, f64::min);
    let max = data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let class_width = ((max - min) / (class_count - 1) as f64 * 100.0).ceil() / 100.0;
    let class_offset = ((min - class_width / 2.0) * 1000.0).floor() / 1000.0;
    (class_width, class_offset)
}

pub fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + a.abs().max(b.abs()))
}

pub fn assert_near(a: f64, b: f64) {
    assert!(near(a, b), "{a} != {b}");
}

/// Series used by the C residue tests: one inner cycle and one open outer cycle.
pub fn curve_series() -> [f64; 4] {
    [2.0, 5.0, 3.0, 7.0]
}

pub fn assert_rows(got: &[[f64; 2]], expect: &[[f64; 2]]) {
    assert_eq!(got.len(), expect.len());
    for (got, expect) in got.iter().zip(expect) {
        assert_near(got[0], expect[0]);
        assert_near(got[1], expect[1]);
    }
}

pub fn assert_vec(got: &[f64], expect: &[f64]) {
    assert_eq!(got.len(), expect.len(), "len {} != {}", got.len(), expect.len());
    for (got, expect) in got.iter().zip(expect) {
        assert_near(*got, *expect);
    }
}

pub fn assert_matrix(got: &Matrix, expect: &Matrix) {
    assert_eq!(got.class_count, expect.class_count);
    assert_vec(&got.counts, &expect.counts);
}

pub fn assert_tp(got: &[TurningPoint], expect: &[TurningPoint]) {
    assert_eq!(got.len(), expect.len());
    for (got, expect) in got.iter().zip(expect) {
        assert_eq!(got.pos, expect.pos);
        assert_eq!(got.adj_pos, expect.adj_pos);
        assert_near(got.value, expect.value);
        assert_near(got.damage, expect.damage);
    }
}
