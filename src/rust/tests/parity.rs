mod common;

use common::{assert_matrix, assert_near, assert_rows, assert_tp, assert_vec, class_param};
use rfcnt::{rfc, LcMethod, Options, ResidualMethod, Rfc};

fn autoresize_matches(use_hcm: bool, use_astm: bool) {
    let data = [1.0, 3.0, 0.5, 10.0, 2.0, 8.0, 1.0];
    let (width, offset) = class_param(&data[..3], 4);
    let mut opt = Options::for_counter();
    opt.class_count = 4;
    opt.class_offset = offset;
    opt.hysteresis = Some(width * 0.99);
    opt.enforce_margin = false;
    opt.residual_method = ResidualMethod::Repeated;
    opt.auto_resize = true;
    opt.use_hcm = use_hcm;
    opt.use_astm = use_astm;

    let one = rfc(&data, width, &opt).unwrap();
    let mut rf = Rfc::new(width, &opt).unwrap();
    rf.feed(&data[..3]).unwrap();
    rf.feed(&data[3..]).unwrap();
    let preview = rf.damage_as(ResidualMethod::Repeated).unwrap();
    let raw_before = rf.res_raw().unwrap();
    rf.finalize(ResidualMethod::Repeated).unwrap();

    assert_near(preview, rf.damage().unwrap());
    assert_near(rf.damage().unwrap(), one.damage);
    assert_vec(&rf.res_raw().unwrap(), &raw_before);
    assert_vec(&rf.res_raw().unwrap(), &one.res_raw);
    assert_rows(&rf.rp().unwrap(), &one.rp);
    assert_matrix(&rf.rfm().unwrap(), &one.rfm);
    assert_tp(&rf.tp().unwrap(), &one.tp);
}

#[test]
fn autoresize_four_point_matches_one_shot() {
    autoresize_matches(false, false);
}

#[test]
fn autoresize_hcm_matches_one_shot() {
    autoresize_matches(true, false);
}

#[test]
fn autoresize_astm_matches_one_shot() {
    autoresize_matches(false, true);
}

#[test]
fn residue_is_readable_before_finalize() {
    let data = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0];
    let (width, offset) = class_param(&data, 4);
    let mut opt = Options::for_counter();
    opt.class_count = 4;
    opt.class_offset = offset;
    opt.hysteresis = Some(width * 0.99);

    let mut rf = Rfc::new(width, &opt).unwrap();
    rf.feed(&data[..3]).unwrap();
    assert!(!rf.residue().unwrap().is_empty());
    let damage_before = rf.damage().unwrap();
    rf.feed(&data[3..]).unwrap();
    assert!(rf.damage().unwrap() + 1e-15 >= damage_before);
    rf.finalize(ResidualMethod::Repeated).unwrap();
}

#[test]
fn close_rejects_further_use() {
    let data = [0.0, 1.0, 0.0, 2.0];
    let (width, offset) = class_param(&data, 10);
    let mut opt = Options::for_counter();
    opt.class_count = 10;
    opt.class_offset = offset;
    opt.hysteresis = Some(width * 0.99);

    let mut rf = Rfc::new(width, &opt).unwrap();
    rf.feed(&data).unwrap();
    rf.close().unwrap();
    assert!(rf.damage().is_err());
    assert!(rf.feed(&[1.0]).is_err());
    assert!(rf.damage_as(ResidualMethod::Repeated).is_err());
}

#[test]
fn fva_t2_counts() {
    let data = [1.0, -2.0, 3.0];
    let expect = [1.0, 1.0, 1.0, 1.0, 1.0, 0.0];
    let mut opt = Options::for_counter();
    opt.class_count = 6;
    opt.class_offset = -2.5;
    opt.hysteresis = Some(0.5);
    opt.enforce_margin = true;
    opt.residual_method = ResidualMethod::Repeated;
    opt.lc_method = LcMethod::Fva;

    let one = rfc(&data, 1.0, &opt).unwrap();
    let counts: Vec<f64> = one.lc.iter().map(|row| row[1]).collect();
    assert_vec(&counts, &expect);

    let mut rf = Rfc::new(1.0, &opt).unwrap();
    rf.feed(&data).unwrap();
    let live: Vec<f64> = rf.lc().unwrap().iter().map(|row| row[1]).collect();
    let preview: Vec<f64> = rf
        .lc_as(ResidualMethod::None)
        .unwrap()
        .iter()
        .map(|row| row[1])
        .collect();
    assert_vec(&live, &expect);
    assert_vec(&preview, &expect);
    rf.finalize(ResidualMethod::Repeated).unwrap();
    let done: Vec<f64> = rf.lc().unwrap().iter().map(|row| row[1]).collect();
    assert_vec(&done, &expect);
}

#[test]
fn fva_zero_class_bound_follows_rising_slope() {
    let data = [-2.0, 2.0, -2.0];
    let mut base = Options::for_counter();
    base.class_count = 5;
    base.class_offset = -2.0;
    base.hysteresis = Some(0.5);
    base.enforce_margin = true;
    base.residual_method = ResidualMethod::Repeated;

    let mut up = base.clone();
    up.lc_method = LcMethod::SlopesUp;
    let mut down = base.clone();
    down.lc_method = LcMethod::SlopesDown;
    let mut fva = base.clone();
    fva.lc_method = LcMethod::Fva;

    let up = rfc(&data, 1.0, &up).unwrap();
    let down = rfc(&data, 1.0, &down).unwrap();
    let fva = rfc(&data, 1.0, &fva).unwrap();

    let mut saw_zero = false;
    for (i, row) in fva.lc.iter().enumerate() {
        let expect = if row[0] >= 0.0 { up.lc[i][1] } else { down.lc[i][1] };
        assert_near(row[1], expect);
        if row[0].abs() <= 1e-12 {
            saw_zero = true;
            assert_near(row[1], up.lc[i][1]);
        }
    }
    assert!(saw_zero);
}
