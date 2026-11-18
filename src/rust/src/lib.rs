//! Rainflow counting for fatigue analysis.
//!
//! The safe API follows the Python package `rfcnt`. Counting is performed by
//! the C core through [`rainflow_ffi.h`](../../lib/rainflow_ffi.h); this crate
//! never sees the packed context struct.

mod ffi;

use std::ffi::CStr;
use std::ptr;

use ffi::{FfiResidue, FfiTp, FfiWl, RfcFfi};

/// Failure from argument checks or from the C core.
#[derive(Debug)]
pub struct Error {
    pub code: i32,
    pub message: String,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for Error {}

fn error(code: i32, message: impl Into<String>) -> Error {
    Error { code, message: message.into() }
}

fn check(code: i32, what: &str) -> Result<(), Error> {
    if code == ffi::OK {
        return Ok(());
    }
    let detail = unsafe { CStr::from_ptr(ffi::rfc_ffi_error_string(code)) }.to_string_lossy();
    Err(error(code, format!("{what} ({detail})")))
}

/// How unclosed cycles are treated by [`rfc`] and [`Rfc::finalize`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum ResidualMethod {
    None = 0,
    Ignore = 1,
    NoFinalize = 2,
    Discard = 3,
    Halfcycles = 4,
    Fullcycles = 5,
    ClormannSeeger = 6,
    Repeated = 7,
    Din45667 = 8,
}

/// How damage increments are spread over the input series. One-shot [`rfc`] only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum SdMethod {
    None = -1,
    Half23 = 0,
    RampAmplitude23 = 1,
    RampDamage23 = 2,
    RampAmplitude24 = 3,
    RampDamage24 = 4,
    FullP2 = 5,
    FullP3 = 6,
    Transient23 = 7,
    Transient23c = 8,
}

/// Which slopes contribute to the level-crossing histogram.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum LcMethod {
    SlopesUp = 0,
    SlopesDown = 1,
    SlopesAll = 2,
    /// FVA Merkblatt. The historical name DIN 45667 is this same value.
    Fva = 3,
}

/// Miner variant used by [`damage_from_rp`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub enum RpDamageCalcMethod {
    Default = 0,
    MinerElementar = 1,
    MinerModified = 2,
    MinerConsistent = 3,
}

/// Wöhler curve passed to [`rfc`] and [`damage_from_rp`].
///
/// [`Woehler::Modified`] matches Python's default `wl_init_modified` path.
/// [`Woehler::Full`] matches `wl_init_any`.
#[derive(Clone, Debug)]
pub enum Woehler {
    Modified { sx: f64, nx: f64, k: f64, k2: f64 },
    Full {
        sx: f64,
        nx: f64,
        k: f64,
        k2: f64,
        sd: f64,
        nd: f64,
        omission: f64,
    },
}

impl Default for Woehler {
    fn default() -> Self {
        Woehler::Modified { sx: 1e3, nx: 1e7, k: 5.0, k2: 5.0 }
    }
}

impl Woehler {
    /// Curve used when [`damage_from_rp`] is called without an explicit curve.
    pub fn for_damage_from_rp() -> Self {
        Woehler::Full {
            sx: 1e3,
            nx: 1e7,
            k: 5.0,
            k2: 5.0,
            sd: 0.0,
            nd: f64::MAX,
            omission: 0.0,
        }
    }
}

/// Miner-consistent (impaired) Wöhler parameters, as returned by Python.
#[derive(Clone, Debug)]
pub struct WoehlerConsistent {
    pub sx: f64,
    pub nx: f64,
    pub sd: f64,
    pub nd: f64,
    pub k: f64,
    pub k2: f64,
    pub q: f64,
    pub q2: f64,
    pub omission: f64,
    pub d: f64,
}

/// Parameters for [`rfc`] and [`Rfc::new`].
///
/// `hysteresis: None` uses `class_width`. [`rfc`] defaults `spread_damage` to
/// [`SdMethod::Transient23c`]. [`Rfc::new`] rejects any other value; use
/// [`Options::for_counter`].
#[derive(Clone, Debug)]
pub struct Options {
    pub class_count: u32,
    pub class_offset: f64,
    pub hysteresis: Option<f64>,
    pub residual_method: ResidualMethod,
    pub spread_damage: SdMethod,
    pub lc_method: LcMethod,
    pub use_hcm: bool,
    pub use_astm: bool,
    pub enforce_margin: bool,
    pub auto_resize: bool,
    pub wl: Woehler,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            class_count: 100,
            class_offset: 0.0,
            hysteresis: None,
            residual_method: ResidualMethod::Repeated,
            spread_damage: SdMethod::Transient23c,
            lc_method: LcMethod::SlopesAll,
            use_hcm: false,
            use_astm: false,
            enforce_margin: true,
            auto_resize: false,
            wl: Woehler::default(),
        }
    }
}

impl Options {
    /// Defaults for [`Rfc`]: damage history stays off.
    pub fn for_counter() -> Self {
        Self {
            spread_damage: SdMethod::None,
            ..Self::default()
        }
    }
}

/// One turning point: position, load, damage, adjacent position. Positions are 1-based.
#[derive(Clone, Debug)]
pub struct TurningPoint {
    pub pos: u64,
    pub value: f64,
    pub damage: f64,
    pub adj_pos: u64,
}

/// Row-major rainflow matrix in full-cycle units (`from` is the row).
#[derive(Clone, Debug)]
pub struct Matrix {
    pub class_count: usize,
    pub counts: Vec<f64>,
}

impl Matrix {
    pub fn get(&self, from: usize, to: usize) -> f64 {
        self.counts[from * self.class_count + to]
    }

    pub fn sum(&self) -> f64 {
        self.counts.iter().sum()
    }
}

/// Result of [`rfc`], with the same fields as the Python dictionary.
#[derive(Clone, Debug)]
pub struct RfcResult {
    pub damage: f64,
    /// Columns are range (`2 * amplitude`) and counts.
    pub rp: Vec<[f64; 2]>,
    /// Columns are class upper bound and counts.
    pub lc: Vec<[f64; 2]>,
    pub tp: Vec<TurningPoint>,
    pub res_raw: Vec<f64>,
    pub res: Vec<f64>,
    pub rfm: Matrix,
    pub dh: Vec<f64>,
    pub wl_miner_consistent: WoehlerConsistent,
}

/// Haigh / FKM amplitude transformation, matching `rfcnt.at_transform`.
#[derive(Clone, Debug)]
pub struct AtOptions {
    pub r_rig: f64,
    pub sm_rig: f64,
    pub r_pinned: bool,
    pub symmetric: bool,
    pub sa_ref: Option<Vec<f64>>,
    pub sm_ref: Option<Vec<f64>>,
}

impl Default for AtOptions {
    fn default() -> Self {
        Self {
            r_rig: -1.0,
            sm_rig: 0.0,
            r_pinned: true,
            symmetric: false,
            sa_ref: None,
            sm_ref: None,
        }
    }
}

struct Session {
    handle: *mut RfcFfi,
}

unsafe impl Send for Session {}

impl Drop for Session {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { ffi::rfc_ffi_destroy(self.handle) };
            self.handle = ptr::null_mut();
        }
    }
}

impl Session {
    fn create() -> Result<Self, Error> {
        let handle = unsafe { ffi::rfc_ffi_create() };
        if handle.is_null() {
            return Err(error(ffi::ERR_INVARG, "Error on memory allocation"));
        }
        Ok(Self { handle })
    }

    fn ptr(&self) -> *mut RfcFfi {
        self.handle
    }

    fn init(&self, class_count: u32, class_width: f64, class_offset: f64, hysteresis: f64) -> Result<(), Error> {
        check(
            unsafe {
                ffi::rfc_ffi_init(
                    self.ptr(),
                    class_count,
                    class_width,
                    class_offset,
                    hysteresis,
                    ffi::FLAGS_DEFAULT,
                )
            },
            "Rainflow initialization error",
        )
    }

    fn apply_flags(&self, opt: &Options) -> Result<(), Error> {
        check(
            unsafe { ffi::rfc_ffi_set_lc_method(self.ptr(), opt.lc_method as i32) },
            "Parameter 'lc_method' must be 0, 1, 2 or 3 (LCMethod enumeration, not a flag mask)!",
        )?;
        let mut flags = 0;
        check(
            unsafe { ffi::rfc_ffi_get_flags(self.ptr(), &mut flags) },
            "Rainflow initialization error",
        )?;
        if opt.auto_resize {
            flags |= ffi::FLAG_AUTORESIZE;
        } else {
            flags &= !ffi::FLAG_AUTORESIZE;
        }
        if opt.enforce_margin {
            flags |= ffi::FLAG_ENFORCE_MARGIN;
        } else {
            flags &= !ffi::FLAG_ENFORCE_MARGIN;
        }
        check(
            unsafe { ffi::rfc_ffi_set_flags(self.ptr(), flags) },
            "Rainflow initialization error",
        )
    }

    fn init_wl(&self, wl: &Woehler) -> Result<(), Error> {
        let code = match *wl {
            Woehler::Modified { sx, nx, k, k2 } => unsafe {
                ffi::rfc_ffi_wl_init_modified(self.ptr(), sx, nx, k, k2)
            },
            Woehler::Full { sx, nx, k, k2, sd, nd, omission } => {
                let param = FfiWl {
                    sd,
                    nd,
                    k: k.abs(),
                    sx,
                    nx,
                    k2: k2.abs(),
                    omission,
                    q: 0.0,
                    q2: 0.0,
                    d: 0.0,
                };
                unsafe { ffi::rfc_ffi_wl_init_any(self.ptr(), &param) }
            }
        };
        check(code, "Rainflow initialization error")
    }

    fn set_counting(&self, opt: &Options) -> Result<(), Error> {
        if opt.use_hcm && opt.use_astm {
            return Err(error(ffi::ERR_INVARG, "`use_HCM` and `use_ASTM` are mutually exclusive!"));
        }
        let method = if opt.use_hcm {
            ffi::COUNT_HCM
        } else if opt.use_astm {
            ffi::COUNT_ASTM
        } else {
            ffi::COUNT_4PTM
        };
        check(
            unsafe { ffi::rfc_ffi_set_counting_method(self.ptr(), method) },
            "Rainflow initialization error",
        )
    }

    fn dh_init(&self, method: SdMethod, cap: u64) -> Result<(), Error> {
        check(
            unsafe { ffi::rfc_ffi_dh_init(self.ptr(), method as i32, cap) },
            "Error on memory allocation",
        )
    }

    fn feed(&self, data: &[f64]) -> Result<(), Error> {
        let ptr = if data.is_empty() { ptr::null() } else { data.as_ptr() };
        check(
            unsafe { ffi::rfc_ffi_feed(self.ptr(), ptr, data.len() as u64) },
            "Error while counting",
        )
    }

    fn finalize(&self, method: ResidualMethod) -> Result<(), Error> {
        check(
            unsafe { ffi::rfc_ffi_finalize(self.ptr(), method as i32) },
            "Error while counting",
        )
    }

    fn class_count(&self) -> Result<u32, Error> {
        let mut n = 0;
        check(
            unsafe { ffi::rfc_ffi_class_count(self.ptr(), &mut n) },
            "Preparing range pair counting",
        )?;
        Ok(n)
    }

    fn damage(&self) -> Result<f64, Error> {
        let mut value = 0.0;
        check(
            unsafe { ffi::rfc_ffi_damage(self.ptr(), &mut value, ptr::null_mut()) },
            "Error reading damage",
        )?;
        Ok(value)
    }

    fn rp_lc(&self) -> Result<(Vec<[f64; 2]>, Vec<[f64; 2]>), Error> {
        let n = self.class_count()? as usize;
        let mut counts = vec![0.0; n];
        let mut axis = vec![0.0; n];
        if n > 0 {
            check(
                unsafe { ffi::rfc_ffi_rp(self.ptr(), counts.as_mut_ptr(), axis.as_mut_ptr()) },
                "Preparing range pair counting",
            )?;
        }
        let rp = counts.iter().zip(axis.iter()).map(|(&c, &sa)| [sa * 2.0, c]).collect();
        if n > 0 {
            check(
                unsafe { ffi::rfc_ffi_lc(self.ptr(), counts.as_mut_ptr(), axis.as_mut_ptr()) },
                "Preparing range pair counting",
            )?;
        }
        let lc = axis.iter().zip(counts.iter()).map(|(&level, &c)| [level, c]).collect();
        Ok((rp, lc))
    }

    fn rfm(&self) -> Result<Matrix, Error> {
        let n = self.class_count()? as usize;
        let mut counts = vec![0.0; n * n];
        if n > 0 {
            check(
                unsafe { ffi::rfc_ffi_rfm(self.ptr(), counts.as_mut_ptr()) },
                "Preparing range pair counting",
            )?;
            let inc = unsafe { ffi::rfc_ffi_full_cycle_increment() };
            if inc != 0.0 {
                for value in &mut counts {
                    *value /= inc;
                }
            }
        }
        Ok(Matrix { class_count: n, counts })
    }

    fn residue_tuples(&self) -> Result<Vec<FfiResidue>, Error> {
        let mut n = 0u32;
        check(
            unsafe { ffi::rfc_ffi_residue(self.ptr(), ptr::null_mut(), &mut n) },
            "Error reading residue",
        )?;
        let mut buf = vec![FfiResidue { value: 0.0, cls: 0 }; n as usize];
        if n > 0 {
            check(
                unsafe { ffi::rfc_ffi_residue(self.ptr(), buf.as_mut_ptr(), &mut n) },
                "Error reading residue",
            )?;
            buf.truncate(n as usize);
        }
        Ok(buf)
    }

    fn turning_points(&self) -> Result<Vec<TurningPoint>, Error> {
        let mut n = 0u64;
        check(
            unsafe { ffi::rfc_ffi_tp(self.ptr(), ptr::null_mut(), &mut n) },
            "Error reading turning points",
        )?;
        let mut buf = vec![
            FfiTp { value: 0.0, damage: 0.0, pos: 0, adj_pos: 0 };
            n as usize
        ];
        if n > 0 {
            check(
                unsafe { ffi::rfc_ffi_tp(self.ptr(), buf.as_mut_ptr(), &mut n) },
                "Error reading turning points",
            )?;
            buf.truncate(n as usize);
        }
        Ok(buf
            .into_iter()
            .map(|tp| TurningPoint {
                pos: tp.pos,
                value: tp.value,
                damage: tp.damage,
                adj_pos: tp.adj_pos,
            })
            .collect())
    }

    fn damage_history(&self, len: usize) -> Result<Vec<f64>, Error> {
        let mut out = vec![0.0; len];
        let mut n = 0u64;
        let code = unsafe { ffi::rfc_ffi_dh(self.ptr(), ptr::null_mut(), &mut n) };
        if code == ffi::ERR_UNSUPPORTED {
            return Ok(out);
        }
        check(code, "Preparing range pair counting")?;
        if n == 0 {
            return Ok(out);
        }
        let mut buf = vec![0.0; n as usize];
        check(
            unsafe { ffi::rfc_ffi_dh(self.ptr(), buf.as_mut_ptr(), &mut n) },
            "Preparing range pair counting",
        )?;
        for (i, value) in buf.into_iter().enumerate() {
            if i < out.len() {
                out[i] = value;
            }
        }
        Ok(out)
    }

    fn wl_consistent(&self) -> Result<WoehlerConsistent, Error> {
        let mut wl = FfiWl {
            sd: 0.0,
            nd: 0.0,
            k: 0.0,
            sx: 0.0,
            nx: 0.0,
            k2: 0.0,
            omission: 0.0,
            q: 0.0,
            q2: 0.0,
            d: 0.0,
        };
        check(
            unsafe { ffi::rfc_ffi_wl_param_get_impaired(self.ptr(), &mut wl) },
            "Error reading Woehler parameters",
        )?;
        Ok(WoehlerConsistent {
            sx: wl.sx,
            nx: wl.nx,
            sd: wl.sd,
            nd: wl.nd,
            k: wl.k,
            k2: wl.k2,
            q: wl.q,
            q2: wl.q2,
            omission: wl.omission,
            d: wl.d,
        })
    }

    fn state(&self) -> i32 {
        unsafe { ffi::rfc_ffi_state(self.ptr()) }
    }

    fn clone_session(&self) -> Result<Session, Error> {
        let mut out = ptr::null_mut();
        check(
            unsafe { ffi::rfc_ffi_clone(self.ptr(), &mut out) },
            "Failed to clone rainflow context",
        )?;
        if out.is_null() {
            return Err(error(ffi::ERR_INVARG, "Failed to clone rainflow context"));
        }
        Ok(Session { handle: out })
    }

    fn at_init(&self, sa_ref: &[f64], sm_ref: &[f64], m: f64, opt: &AtOptions) -> Result<(), Error> {
        if m < 0.0 {
            return Err(error(ffi::ERR_INVARG, "`M` (mean stress sensitivity) must be >= 0"));
        }
        let has_ref = opt.sa_ref.is_some() || opt.sm_ref.is_some();
        if has_ref && (opt.sa_ref.is_none() || opt.sm_ref.is_none()) {
            return Err(error(
                ffi::ERR_INVARG,
                "`Sa_ref` and `Sm_ref` must both be given or both omitted",
            ));
        }
        if has_ref && (sa_ref.len() != sm_ref.len() || sa_ref.len() < 2) {
            return Err(error(
                ffi::ERR_INVARG,
                "`Sa_ref` and `Sm_ref` must be 1-D arrays of equal length",
            ));
        }
        if has_ref && opt.symmetric {
            return Err(error(
                ffi::ERR_INVARG,
                "`symmetric` is not supported with a custom Haigh reference curve",
            ));
        }
        let (sa_ptr, sm_ptr, count) = if sa_ref.is_empty() {
            (ptr::null(), ptr::null(), 0)
        } else {
            (sa_ref.as_ptr(), sm_ref.as_ptr(), sa_ref.len() as u32)
        };
        check(
            unsafe {
                ffi::rfc_ffi_at_init(
                    self.ptr(),
                    sa_ptr,
                    sm_ptr,
                    count,
                    m,
                    opt.sm_rig,
                    opt.r_rig,
                    i32::from(opt.r_pinned),
                    i32::from(opt.symmetric),
                )
            },
            "Amplitude transformation initialization error",
        )
    }

    fn at_transform(&self, sa: &[f64], sm: &[f64]) -> Result<Vec<f64>, Error> {
        if sa.len() != sm.len() {
            return Err(error(ffi::ERR_INVARG, "`Sa` and `Sm` must have the same number of elements"));
        }
        let mut out = Vec::with_capacity(sa.len());
        for (&sa_i, &sm_i) in sa.iter().zip(sm.iter()) {
            let mut sa_t = 0.0;
            check(
                unsafe { ffi::rfc_ffi_at_transform(self.ptr(), sa_i, sm_i, &mut sa_t) },
                "Error while amplitude transformation",
            )?;
            out.push(sa_t);
        }
        Ok(out)
    }
}

fn open_configured(class_width: f64, opt: &Options) -> Result<Session, Error> {
    let hysteresis = opt.hysteresis.unwrap_or(class_width);
    let session = Session::create()?;
    session.init(opt.class_count, class_width, opt.class_offset, hysteresis)?;
    session.apply_flags(opt)?;
    session.init_wl(&opt.wl)?;
    session.set_counting(opt)?;
    Ok(session)
}

/// Drop a pending closed 4-point cycle so `res_raw` contains only open cycles.
/// Same rule as `strip_pending_closed_cycle` in the Python extension.
fn strip_pending_closed_cycle(residue: &mut Vec<FfiResidue>) {
    if residue.len() < 4 {
        return;
    }
    let idx = residue.len() - 4;
    let mut a = residue[idx].cls;
    let mut b = residue[idx + 1].cls;
    let mut c = residue[idx + 2].cls;
    let mut d = residue[idx + 3].cls;
    if b > c {
        std::mem::swap(&mut b, &mut c);
    }
    if a > d {
        std::mem::swap(&mut a, &mut d);
    }
    if a <= b && c <= d {
        let len = residue.len();
        residue.drain(len - 3..len - 1);
    }
}

fn res_raw_from(session: &Session) -> Result<Vec<f64>, Error> {
    let mut tuples = session.residue_tuples()?;
    strip_pending_closed_cycle(&mut tuples);
    Ok(tuples.into_iter().map(|point| point.value).collect())
}

fn snapshot(session: &Session, dh_len: usize) -> Result<RfcResult, Error> {
    let damage = session.damage()?;
    let (rp, lc) = session.rp_lc()?;
    let tp = session.turning_points()?;
    let res_raw = res_raw_from(session)?;
    let res = session.residue_tuples()?.into_iter().map(|point| point.value).collect();
    let rfm = session.rfm()?;
    let dh = session.damage_history(dh_len)?;
    let wl_miner_consistent = session.wl_consistent()?;
    Ok(RfcResult {
        damage,
        rp,
        lc,
        tp,
        res_raw,
        res,
        rfm,
        dh,
        wl_miner_consistent,
    })
}

/// Count `data` in one shot. Defaults match `rfcnt.rfc`.
pub fn rfc(data: &[f64], class_width: f64, opt: &Options) -> Result<RfcResult, Error> {
    let session = open_configured(class_width, opt)?;
    if opt.spread_damage != SdMethod::None {
        session.dh_init(opt.spread_damage, data.len() as u64)?;
    }
    session.feed(data)?;
    let res_raw = res_raw_from(&session)?;
    session.finalize(opt.residual_method)?;
    let mut result = snapshot(&session, data.len())?;
    result.res_raw = res_raw;
    Ok(result)
}

/// Stateful counter. Chunked [`Rfc::feed`] cannot keep a damage-history stream;
/// use [`rfc`] when `spread_damage` is not [`SdMethod::None`].
pub struct Rfc {
    session: Option<Session>,
    finalized: bool,
    res_raw: Option<Vec<f64>>,
    at_sa: Vec<f64>,
    at_sm: Vec<f64>,
}

impl Rfc {
    pub fn new(class_width: f64, opt: &Options) -> Result<Self, Error> {
        if opt.spread_damage != SdMethod::None {
            return Err(error(
                ffi::ERR_INVARG,
                "Damage history (spread_damage) is not supported on the stateful RFC class; \
                 chunked feed() cannot keep the original sample stream. Use rfc() for \
                 one-shot counting with damage history.",
            ));
        }
        let session = open_configured(class_width, opt)?;
        Ok(Self {
            session: Some(session),
            finalized: false,
            res_raw: None,
            at_sa: Vec::new(),
            at_sm: Vec::new(),
        })
    }

    fn ensure_feedable(&self) -> Result<&Session, Error> {
        if self.finalized || self.session.is_none() {
            return Err(error(ffi::ERR_INVARG, "RFC object is finalized or closed"));
        }
        Ok(self.session.as_ref().unwrap())
    }

    fn ensure_alive(&self) -> Result<&Session, Error> {
        match &self.session {
            Some(session) if session.state() >= ffi::STATE_INIT => Ok(session),
            _ => Err(error(ffi::ERR_INVARG, "RFC object is closed")),
        }
    }

    pub fn feed(&mut self, data: &[f64]) -> Result<(), Error> {
        self.ensure_feedable()?.feed(data)
    }

    pub fn finalize(&mut self, residual_method: ResidualMethod) -> Result<(), Error> {
        let session = self.ensure_feedable()?;
        let raw = res_raw_from(session)?;
        if let Err(err) = session.finalize(residual_method) {
            return Err(error(err.code, err.message.replace("Error while counting", "Error while finalizing")));
        }
        self.res_raw = Some(raw);
        self.finalized = true;
        Ok(())
    }

    pub fn close(&mut self) -> Result<(), Error> {
        if let Some(session) = self.session.as_ref() {
            if session.state() >= ffi::STATE_INIT {
                check(unsafe { ffi::rfc_ffi_deinit(session.ptr()) }, "Error while closing")?;
            }
        }
        self.finalized = true;
        Ok(())
    }

    pub fn state(&self) -> Result<i32, Error> {
        Ok(self.ensure_alive()?.state())
    }

    pub fn damage(&self) -> Result<f64, Error> {
        self.ensure_alive()?.damage()
    }

    pub fn residue(&self) -> Result<Vec<f64>, Error> {
        Ok(self.ensure_alive()?.residue_tuples()?.into_iter().map(|p| p.value).collect())
    }

    pub fn rp(&self) -> Result<Vec<[f64; 2]>, Error> {
        Ok(self.ensure_alive()?.rp_lc()?.0)
    }

    pub fn lc(&self) -> Result<Vec<[f64; 2]>, Error> {
        Ok(self.ensure_alive()?.rp_lc()?.1)
    }

    pub fn rfm(&self) -> Result<Matrix, Error> {
        self.ensure_alive()?.rfm()
    }

    pub fn tp(&self) -> Result<Vec<TurningPoint>, Error> {
        self.ensure_alive()?.turning_points()
    }

    pub fn res_raw(&self) -> Result<Vec<f64>, Error> {
        if let Some(frozen) = &self.res_raw {
            return Ok(frozen.clone());
        }
        res_raw_from(self.ensure_alive()?)
    }

    pub fn wl_miner_consistent(&self) -> Result<WoehlerConsistent, Error> {
        self.ensure_alive()?.wl_consistent()
    }

    fn preview(&self, residual_method: ResidualMethod) -> Result<Session, Error> {
        let session = self.ensure_feedable()?;
        let clone = session.clone_session()?;
        clone.finalize(residual_method).map_err(|err| {
            error(err.code, err.message.replace("Error while counting", "Error while preview finalize"))
        })?;
        Ok(clone)
    }

    pub fn damage_as(&self, residual_method: ResidualMethod) -> Result<f64, Error> {
        self.preview(residual_method)?.damage()
    }

    pub fn rp_as(&self, residual_method: ResidualMethod) -> Result<Vec<[f64; 2]>, Error> {
        Ok(self.preview(residual_method)?.rp_lc()?.0)
    }

    pub fn lc_as(&self, residual_method: ResidualMethod) -> Result<Vec<[f64; 2]>, Error> {
        Ok(self.preview(residual_method)?.rp_lc()?.1)
    }

    pub fn rfm_as(&self, residual_method: ResidualMethod) -> Result<Matrix, Error> {
        self.preview(residual_method)?.rfm()
    }

    /// Install the Haigh transformation. Call this after [`Rfc::new`] and before [`Rfc::feed`].
    ///
    /// Reference-curve buffers are kept alive on this counter. A preview clone aliases them.
    pub fn at_init(&mut self, m: f64, opt: &AtOptions) -> Result<(), Error> {
        let session = self.ensure_feedable()?;
        if session.state() != ffi::STATE_INIT {
            return Err(error(
                ffi::ERR_INVARG,
                "at_init() must be called after construction and before feed()",
            ));
        }
        let sa_ref = opt.sa_ref.clone().unwrap_or_default();
        let sm_ref = opt.sm_ref.clone().unwrap_or_default();
        // The core stores these pointers. Moving the Vec afterwards keeps the allocation.
        session.at_init(&sa_ref, &sm_ref, m, opt)?;
        self.at_sa = sa_ref;
        self.at_sm = sm_ref;
        Ok(())
    }

    pub fn at_transform(&self, sa: &[f64], sm: &[f64]) -> Result<Vec<f64>, Error> {
        self.ensure_alive()?.at_transform(sa, sm)
    }
}

/// Damage from a range-pair histogram. `Sa` is sorted ascending, matching Python.
///
/// `wl: None` uses the Python default curve (`sx=1e3`, `nx=1e7`, `k=k2=5`).
pub fn damage_from_rp(
    sa: &[f64],
    counts: &[f64],
    wl: Option<&Woehler>,
    method: RpDamageCalcMethod,
) -> Result<f64, Error> {
    if (method as i32) < 0 || method as i32 > 3 {
        return Err(error(ffi::ERR_INVARG, "`method` must be in range 0 to 3."));
    }
    if sa.len() != counts.len() {
        return Err(error(ffi::ERR_INVARG, "`Sa` and `counts` must be 1-D arrays of equal length"));
    }
    if counts.iter().any(|count| *count < 0.0) {
        return Err(error(ffi::ERR_INVARG, "Negative values in `counts`."));
    }
    let mut order: Vec<usize> = (0..sa.len()).collect();
    order.sort_by(|&a, &b| sa[a].partial_cmp(&sa[b]).unwrap_or(std::cmp::Ordering::Equal));
    let sa_sorted: Vec<f64> = order.iter().map(|&i| sa[i]).collect();
    let counts_sorted: Vec<f64> = order.iter().map(|&i| counts[i]).collect();
    let curve = wl.cloned().unwrap_or_else(Woehler::for_damage_from_rp);
    let session = Session::create()?;
    session.init(sa_sorted.len() as u32, 1.0, 0.0, 0.0)?;
    session.init_wl(&curve)?;
    let mut damage = 0.0;
    check(
        unsafe {
            ffi::rfc_ffi_damage_from_rp(
                session.ptr(),
                &mut damage,
                counts_sorted.as_ptr(),
                sa_sorted.as_ptr(),
                method as i32,
            )
        },
        "Error while calculation",
    )?;
    Ok(damage)
}

/// One-shot Haigh transformation. `M` is the mean-stress sensitivity.
pub fn at_transform(sa: &[f64], sm: &[f64], m: f64, opt: &AtOptions) -> Result<Vec<f64>, Error> {
    let session = Session::create()?;
    session.init(0, 0.0, 0.0, 0.0)?;
    let sa_ref = opt.sa_ref.as_deref().unwrap_or(&[]);
    let sm_ref = opt.sm_ref.as_deref().unwrap_or(&[]);
    session.at_init(sa_ref, sm_ref, m, opt)?;
    session.at_transform(sa, sm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn class_param(data: &[f64], class_count: u32) -> (f64, f64) {
        if data.is_empty() {
            return (1.0, 0.0);
        }
        let min = data.iter().copied().fold(f64::INFINITY, f64::min);
        let max = data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let class_width = ((max - min) / (class_count - 1) as f64 * 100.0).ceil() / 100.0;
        let class_offset = ((min - class_width / 2.0) * 1000.0).floor() / 1000.0;
        (class_width, class_offset)
    }

    fn near(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9
    }

    fn series_options(data: &[f64], class_count: u32, hyst_scale: f64) -> (f64, Options) {
        let (class_width, class_offset) = class_param(data, class_count);
        let mut opt = Options::for_counter();
        opt.class_count = class_count;
        opt.class_offset = class_offset;
        opt.hysteresis = Some(class_width * hyst_scale);
        opt.enforce_margin = false;
        opt.residual_method = ResidualMethod::None;
        (class_width, opt)
    }

    #[test]
    fn empty_series() {
        let (width, opt) = series_options(&[], 100, 1.0);
        let res = rfc(&[], width, &opt).unwrap();
        assert!(near(res.rfm.sum(), 0.0));
        assert!(res.res.is_empty());
    }

    #[test]
    fn single_cycle_up() {
        let data = [1.0, 3.0, 2.0, 4.0];
        let (width, opt) = series_options(&data, 4, 0.99);
        let res = rfc(&data, width, &opt).unwrap();
        assert!(near(res.rfm.sum(), 1.0));
        assert!(near(res.rfm.get(2, 1), 1.0));
        assert_eq!(res.res, vec![1.0, 4.0]);
    }

    #[test]
    fn single_cycle_down() {
        let data = [4.0, 2.0, 3.0, 1.0];
        let (width, opt) = series_options(&data, 4, 0.99);
        let res = rfc(&data, width, &opt).unwrap();
        assert!(near(res.rfm.sum(), 1.0));
        assert!(near(res.rfm.get(1, 2), 1.0));
        assert_eq!(res.res, vec![4.0, 1.0]);
    }

    #[test]
    fn small_sample() {
        let data = [2.0, 5.0, 3.0, 6.0, 2.0, 4.0, 1.0, 6.0, 1.0, 4.0, 1.0, 5.0, 3.0, 6.0, 3.0, 6.0, 1.0, 5.0, 2.0];
        let (width, opt) = series_options(&data, 6, 1.0);
        let res = rfc(&data, width, &opt).unwrap();
        assert!(near(res.rfm.sum(), 7.0));
        assert!(near(res.rfm.get(4, 2), 2.0));
        assert!(near(res.rfm.get(5, 2), 1.0));
        assert!(near(res.rfm.get(0, 3), 1.0));
        assert!(near(res.rfm.get(1, 3), 1.0));
        assert!(near(res.rfm.get(0, 5), 2.0));
        assert_eq!(res.res, vec![2.0, 6.0, 1.0, 5.0, 2.0]);
    }

    #[test]
    fn chunked_feed_matches_one_shot() {
        let data = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0];
        let (width, mut opt) = series_options(&data, 4, 0.99);
        opt.residual_method = ResidualMethod::Repeated;
        let one = rfc(&data, width, &opt).unwrap();
        let mut rf = Rfc::new(width, &opt).unwrap();
        rf.feed(&data[..4]).unwrap();
        rf.feed(&data[4..]).unwrap();
        rf.finalize(ResidualMethod::Repeated).unwrap();
        assert!(near(rf.damage().unwrap(), one.damage));
        assert_eq!(rf.residue().unwrap(), one.res);
        assert_eq!(rf.rp().unwrap(), one.rp);
        assert_eq!(rf.lc().unwrap(), one.lc);
        assert_eq!(rf.rfm().unwrap().counts, one.rfm.counts);
        assert_eq!(rf.res_raw().unwrap(), one.res_raw);
        assert!(near(rf.wl_miner_consistent().unwrap().sx, one.wl_miner_consistent.sx));
        assert!(near(rf.wl_miner_consistent().unwrap().k, one.wl_miner_consistent.k));
    }

    #[test]
    fn preview_matches_finalize_and_does_not_mutate() {
        let data = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0];
        let (width, mut opt) = series_options(&data, 4, 0.99);
        opt.residual_method = ResidualMethod::Repeated;
        let one = rfc(&data, width, &opt).unwrap();
        let mut rf = Rfc::new(width, &opt).unwrap();
        rf.feed(&data).unwrap();
        let residue_before = rf.residue().unwrap();
        let damage_before = rf.damage().unwrap();
        let state_before = rf.state().unwrap();
        let preview = rf.damage_as(ResidualMethod::Repeated).unwrap();
        let rp_preview = rf.rp_as(ResidualMethod::Repeated).unwrap();
        let lc_preview = rf.lc_as(ResidualMethod::Repeated).unwrap();
        let rfm_preview = rf.rfm_as(ResidualMethod::Repeated).unwrap();
        assert_eq!(rf.state().unwrap(), state_before);
        assert_eq!(rf.residue().unwrap(), residue_before);
        assert!(near(rf.damage().unwrap(), damage_before));
        rf.finalize(ResidualMethod::Repeated).unwrap();
        assert!(near(rf.damage().unwrap(), preview));
        assert_eq!(rf.rp().unwrap(), rp_preview);
        assert_eq!(rf.lc().unwrap(), lc_preview);
        assert_eq!(rf.rfm().unwrap().counts, rfm_preview.counts);
        assert!(near(preview, one.damage));
    }

    #[test]
    fn feed_and_preview_after_finalize_fail() {
        let data = [0.0, 1.0, 0.5, 2.0];
        let (width, opt) = series_options(&data, 10, 0.99);
        let mut rf = Rfc::new(width, &opt).unwrap();
        rf.feed(&data).unwrap();
        rf.finalize(ResidualMethod::Repeated).unwrap();
        assert!(rf.feed(&[1.0]).is_err());
        assert!(rf.damage_as(ResidualMethod::Repeated).is_err());
    }

    #[test]
    fn hcm_and_astm_match_one_shot() {
        let data = [1.0, 3.0, 2.0, 4.0, 1.0, 5.0, 2.0, 4.0];
        for (hcm, astm) in [(true, false), (false, true)] {
            let (width, mut opt) = series_options(&data, 4, 0.99);
            opt.residual_method = ResidualMethod::Repeated;
            opt.use_hcm = hcm;
            opt.use_astm = astm;
            let one = rfc(&data, width, &opt).unwrap();
            let mut rf = Rfc::new(width, &opt).unwrap();
            rf.feed(&data[..4]).unwrap();
            rf.feed(&data[4..]).unwrap();
            let preview = rf.damage_as(ResidualMethod::Repeated).unwrap();
            let raw_before = rf.res_raw().unwrap();
            rf.finalize(ResidualMethod::Repeated).unwrap();
            assert_eq!(rf.res_raw().unwrap(), raw_before);
            assert!(near(preview, one.damage));
            assert!(near(rf.damage().unwrap(), one.damage));
            assert_eq!(rf.rfm().unwrap().counts, one.rfm.counts);
            assert_eq!(rf.tp().unwrap().len(), one.tp.len());
        }
    }

    #[test]
    fn hcm_and_astm_are_exclusive() {
        let mut opt = Options::for_counter();
        opt.use_hcm = true;
        opt.use_astm = true;
        assert!(Rfc::new(1.0, &opt).is_err());
    }

    #[test]
    fn spread_damage_rejected_on_counter() {
        let opt = Options::default();
        assert!(Rfc::new(1.0, &opt).is_err());
    }

    #[test]
    fn level_crossing_fva_and_slopes() {
        let data = [-3.0, 2.0, -1.0];
        let mut opt = Options::for_counter();
        opt.class_count = 6;
        opt.class_offset = -3.5;
        opt.hysteresis = Some(0.5);
        opt.enforce_margin = true;
        opt.residual_method = ResidualMethod::Repeated;
        opt.lc_method = LcMethod::Fva;
        let res = rfc(&data, 1.0, &opt).unwrap();
        let levels: Vec<f64> = res.lc.iter().map(|row| row[0]).collect();
        let counts: Vec<f64> = res.lc.iter().map(|row| row[1]).collect();
        assert_eq!(levels, vec![-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]);
        assert_eq!(counts, vec![0.0, 0.0, 1.0, 1.0, 1.0, 0.0]);

        let mut both = opt.clone();
        both.lc_method = LcMethod::SlopesAll;
        let mut up = opt.clone();
        up.lc_method = LcMethod::SlopesUp;
        let x = [0.0, 1.0, 0.0, 2.0];
        let (width, offset) = class_param(&x, 10);
        for item in [&mut both, &mut up] {
            item.class_count = 10;
            item.class_offset = offset;
            item.hysteresis = Some(width * 0.99);
            item.enforce_margin = false;
        }
        let all = rfc(&x, width, &both).unwrap();
        let mut default_opt = both.clone();
        default_opt.lc_method = LcMethod::SlopesAll;
        let default = rfc(&x, width, &default_opt).unwrap();
        let rising = rfc(&x, width, &up).unwrap();
        assert_eq!(default.lc, all.lc);
        assert_ne!(default.lc, rising.lc);
        assert_eq!(LcMethod::Fva as i32, 3);
    }

    #[test]
    fn fva_follows_sign_of_class_bound() {
        let x = [-3.0, 2.0, -1.0, 2.5, -2.0];
        let mut base = Options::for_counter();
        base.class_count = 7;
        base.class_offset = -3.5;
        base.hysteresis = Some(0.5);
        base.enforce_margin = true;
        base.residual_method = ResidualMethod::Repeated;
        let mut up = base.clone();
        up.lc_method = LcMethod::SlopesUp;
        let mut dn = base.clone();
        dn.lc_method = LcMethod::SlopesDown;
        let mut fva = base.clone();
        fva.lc_method = LcMethod::Fva;
        let up = rfc(&x, 1.0, &up).unwrap();
        let dn = rfc(&x, 1.0, &dn).unwrap();
        let fva = rfc(&x, 1.0, &fva).unwrap();
        for (i, row) in fva.lc.iter().enumerate() {
            let expect = if row[0] >= 0.0 { up.lc[i][1] } else { dn.lc[i][1] };
            assert!(near(row[1], expect));
        }
        assert_ne!(up.lc, fva.lc);
    }

    #[test]
    fn at_transform_matches_python_reference() {
        let sa = [0.0, 0.1, 1.0, 2.0, 2.0, 3.0, 3.0, 2.0, 3.0, 4.0, 3.0, 2.0, 2.0, 0.2];
        let sm = [2.0, 9.0, 4.0, 4.0, 2.0, 3.0, 2.0, 1.0, 1.0, 1.0, 0.0, -2.0, -9.0, -9.0];
        let expected = [0.0, 0.153636, 1.536363, 2.836363, 2.6, 3.9, 3.6, 2.3, 3.3, 4.3, 3.0, 1.4, 1.4, 0.14];
        let got = at_transform(&sa, &sm, 0.3, &AtOptions::default()).unwrap();
        for (got, exp) in got.iter().zip(expected.iter()) {
            assert!((got - exp).abs() <= 1e-5, "{got} != {exp}");
        }
        let flat = at_transform(&[3.0], &[0.0], 0.3, &AtOptions::default()).unwrap();
        assert!((flat[0] - 3.0).abs() <= 1e-10);
        let mut sym = AtOptions::default();
        sym.symmetric = true;
        let got_sym = at_transform(&[3.0], &[1.0], 0.3, &sym).unwrap();
        assert!((got_sym[0] - 3.3).abs() <= 1e-10);

        let mut rf = Rfc::new(1.0, &Options::for_counter()).unwrap();
        rf.at_init(0.3, &AtOptions::default()).unwrap();
        let got = rf.at_transform(&[3.0, 2.0], &[1.0, 2.0]).unwrap();
        assert!((got[0] - 3.3).abs() <= 1e-10);
        assert!((got[1] - 2.6).abs() <= 1e-10);
        rf.feed(&[0.0, 1.0, 0.0]).unwrap();
        assert!(rf.at_init(0.3, &AtOptions::default()).is_err());
    }

    #[test]
    fn damage_from_rp_is_positive_and_rejects_bad_input() {
        let damage = damage_from_rp(&[10.0, 100.0], &[1.0, 2.0], None, RpDamageCalcMethod::MinerElementar).unwrap();
        assert!(damage > 0.0);
        assert!(damage_from_rp(&[1.0], &[-1.0], None, RpDamageCalcMethod::Default).is_err());
        assert!(damage_from_rp(&[1.0], &[1.0, 2.0], None, RpDamageCalcMethod::Default).is_err());
    }
}
