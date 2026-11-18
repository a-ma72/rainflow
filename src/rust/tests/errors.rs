use rfcnt::{Options, ResidualMethod, Rfc};

#[test]
fn hcm_and_astm_together_are_rejected() {
    let mut opt = Options::for_counter();
    opt.use_hcm = true;
    opt.use_astm = true;
    let Err(err) = Rfc::new(1.0, &opt) else {
        panic!("HCM and ASTM together must be rejected");
    };
    assert!(err.message.contains("mutually exclusive"), "{}", err.message);
}

#[test]
fn spread_damage_on_counter_mentions_rfc() {
    let Err(err) = Rfc::new(1.0, &Options::default()) else {
        panic!("spread_damage on Rfc must be rejected");
    };
    assert!(err.message.contains("rfc()"), "{}", err.message);
}

#[test]
fn feed_after_finalize_fails() {
    let data = [0.0, 1.0, 0.5, 2.0];
    let mut opt = Options::for_counter();
    opt.class_count = 10;
    let mut rf = Rfc::new(1.0, &opt).unwrap();
    rf.feed(&data).unwrap();
    rf.finalize(ResidualMethod::Repeated).unwrap();
    assert!(rf.feed(&[1.0]).is_err());
}

#[test]
fn feed_after_close_fails() {
    let mut rf = Rfc::new(1.0, &Options::for_counter()).unwrap();
    rf.feed(&[0.0, 1.0, 0.0]).unwrap();
    rf.close().unwrap();
    assert!(rf.feed(&[1.0]).is_err());
}

#[test]
fn second_finalize_fails() {
    let data = [0.0, 1.0, 0.0, 2.0];
    let mut rf = Rfc::new(1.0, &Options::for_counter()).unwrap();
    rf.feed(&data).unwrap();
    rf.finalize(ResidualMethod::Repeated).unwrap();
    let Err(err) = rf.finalize(ResidualMethod::Repeated) else {
        panic!("second finalize must fail");
    };
    assert!(err.message.contains("finalized"), "{}", err.message);
}
