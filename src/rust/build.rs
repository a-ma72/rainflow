fn main() {
    let manifest = std::path::PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let lib = manifest.join("../lib");

    let mut build = cc::Build::new();
    build
        .file(lib.join("rainflow.c"))
        .file(lib.join("rainflow_ffi.c"))
        .include(&lib)
        .std("c11")
        .warnings(false)
        .define("RFC_USE_INTEGRAL_COUNTS", "0")
        .define("RFC_USE_HYSTERESIS_FILTER", "1")
        .define("RFC_MINIMAL", "0")
        .define("RFC_TP_SUPPORT", "1")
        .define("RFC_HCM_SUPPORT", "1")
        .define("RFC_ASTM_SUPPORT", "1")
        .define("RFC_USE_DELEGATES", "1")
        .define("RFC_GLOBAL_EXTREMA", "1")
        .define("RFC_DAMAGE_FAST", "1")
        .define("RFC_DH_SUPPORT", "1")
        .define("RFC_AT_SUPPORT", "1")
        .define("RFC_AR_SUPPORT", "1")
        .define("RFC_DEBUG_FLAGS", "0")
        .define("RFC_EXPORT_MEX", "0")
        .define("RFC_EXPORT_PY", "0");

    if std::env::var("CARGO_CFG_TARGET_ENV").ok().as_deref() == Some("msvc") {
        build.define("_CRT_SECURE_NO_WARNINGS", None);
    }

    build.compile("rainflow");

    println!("cargo:rerun-if-changed=../lib/rainflow.c");
    println!("cargo:rerun-if-changed=../lib/rainflow.h");
    println!("cargo:rerun-if-changed=../lib/rainflow_ffi.c");
    println!("cargo:rerun-if-changed=../lib/rainflow_ffi.h");

    if std::env::var("CARGO_CFG_TARGET_OS").ok().as_deref() != Some("windows") {
        println!("cargo:rustc-link-lib=m");
    }
}
