//! Hand-written bindings for `rainflow_ffi.h`. The packed `rfc_ctx_s` never appears here.

#![allow(dead_code)]

use std::os::raw::c_char;

pub const OK: i32 = 0;
pub const ERR_INVARG: i32 = 1;
pub const ERR_UNSUPPORTED: i32 = 2;

pub const RES_NONE: i32 = 0;
pub const RES_REPEATED: i32 = 7;

pub const SD_NONE: i32 = -1;
pub const SD_TRANSIENT_23C: i32 = 8;

pub const LC_SLOPES_UP: i32 = 0;
pub const LC_SLOPES_DOWN: i32 = 1;
pub const LC_SLOPES_ALL: i32 = 2;
pub const LC_FVA: i32 = 3;

pub const RP_DEFAULT: i32 = 0;

pub const COUNT_4PTM: i32 = 1;
pub const COUNT_HCM: i32 = 2;
pub const COUNT_ASTM: i32 = 3;

pub const FLAG_COUNT_LC: i32 = (1 << 4) | (1 << 5);
pub const FLAG_ENFORCE_MARGIN: i32 = 1 << 7;
pub const FLAG_AUTORESIZE: i32 = 1 << 11;
pub const FLAGS_DEFAULT: i32 = -1;

pub const STATE_INIT: i32 = 1;

#[repr(C)]
pub struct RfcFfi {
    _private: [u8; 0],
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct FfiWl {
    pub sd: f64,
    pub nd: f64,
    pub k: f64,
    pub sx: f64,
    pub nx: f64,
    pub k2: f64,
    pub omission: f64,
    pub q: f64,
    pub q2: f64,
    pub d: f64,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct FfiTp {
    pub value: f64,
    pub damage: f64,
    pub pos: u64,
    pub adj_pos: u64,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct FfiResidue {
    pub value: f64,
    pub cls: u32,
}

extern "C" {
    pub fn rfc_ffi_create() -> *mut RfcFfi;
    pub fn rfc_ffi_destroy(h: *mut RfcFfi);
    pub fn rfc_ffi_deinit(h: *mut RfcFfi) -> i32;
    pub fn rfc_ffi_init(
        h: *mut RfcFfi,
        class_count: u32,
        class_width: f64,
        class_offset: f64,
        hysteresis: f64,
        flags: i32,
    ) -> i32;
    pub fn rfc_ffi_set_flags(h: *mut RfcFfi, flags: i32) -> i32;
    pub fn rfc_ffi_get_flags(h: *const RfcFfi, flags: *mut i32) -> i32;
    pub fn rfc_ffi_set_lc_method(h: *mut RfcFfi, method: i32) -> i32;
    pub fn rfc_ffi_set_counting_method(h: *mut RfcFfi, method: i32) -> i32;
    pub fn rfc_ffi_wl_init_modified(h: *mut RfcFfi, sx: f64, nx: f64, k: f64, k2: f64) -> i32;
    pub fn rfc_ffi_wl_init_any(h: *mut RfcFfi, wl: *const FfiWl) -> i32;
    pub fn rfc_ffi_wl_param_get_impaired(h: *const RfcFfi, out: *mut FfiWl) -> i32;
    pub fn rfc_ffi_dh_init(h: *mut RfcFfi, method: i32, cap: u64) -> i32;
    pub fn rfc_ffi_feed(h: *mut RfcFfi, data: *const f64, count: u64) -> i32;
    pub fn rfc_ffi_finalize(h: *mut RfcFfi, residual_method: i32) -> i32;
    pub fn rfc_ffi_class_count(h: *const RfcFfi, count: *mut u32) -> i32;
    pub fn rfc_ffi_damage(h: *const RfcFfi, damage: *mut f64, damage_residue: *mut f64) -> i32;
    pub fn rfc_ffi_rp(h: *const RfcFfi, counts: *mut f64, amplitude: *mut f64) -> i32;
    pub fn rfc_ffi_lc(h: *const RfcFfi, counts: *mut f64, level: *mut f64) -> i32;
    pub fn rfc_ffi_rfm(h: *const RfcFfi, rfm_row_major: *mut f64) -> i32;
    pub fn rfc_ffi_residue(h: *const RfcFfi, buf: *mut FfiResidue, count: *mut u32) -> i32;
    pub fn rfc_ffi_tp(h: *const RfcFfi, buf: *mut FfiTp, count: *mut u64) -> i32;
    pub fn rfc_ffi_dh(h: *const RfcFfi, buf: *mut f64, count: *mut u64) -> i32;
    pub fn rfc_ffi_state(h: *const RfcFfi) -> i32;
    pub fn rfc_ffi_error_string(err: i32) -> *const c_char;
    pub fn rfc_ffi_full_cycle_increment() -> f64;
    pub fn rfc_ffi_clone(src: *const RfcFfi, out: *mut *mut RfcFfi) -> i32;
    pub fn rfc_ffi_at_init(
        h: *mut RfcFfi,
        sa: *const f64,
        sm: *const f64,
        count: u32,
        m: f64,
        sm_rig: f64,
        r_rig: f64,
        r_pinned: i32,
        symmetric: i32,
    ) -> i32;
    pub fn rfc_ffi_at_transform(h: *const RfcFfi, sa: f64, sm: f64, sa_out: *mut f64) -> i32;
    pub fn rfc_ffi_damage_from_rp(
        h: *mut RfcFfi,
        damage: *mut f64,
        counts: *const f64,
        sa: *const f64,
        method: i32,
    ) -> i32;
}
