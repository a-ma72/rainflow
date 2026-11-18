// rfcnt.cpp — pybind11 port of the rfcnt Python extension.
//
// Argument parsing, keyword-only defaults, GIL handling and the NumPy/dict
// return shapes stay here. Counting goes through rainflow_ffi.h, the same
// opaque handle the Rust crate uses. Turning points, the Wöhler curve and
// the preview clone live inside that handle.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>   // std::optional <-> None conversion

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <rainflow_ffi.h>

#include "docstrings_generated.hpp"

namespace py = pybind11;

// Owns one rfc_ffi handle. close() deinits without freeing; the destructor
// always destroys, and destroy is a no-op on a handle that is already deinited.
struct FfiHandle {
    rfc_ffi* h = nullptr;

    explicit FfiHandle(rfc_ffi* p = nullptr) : h(p) {}
    ~FfiHandle() { rfc_ffi_destroy(h); }

    FfiHandle(const FfiHandle&) = delete;
    FfiHandle& operator=(const FfiHandle&) = delete;

    FfiHandle(FfiHandle&& other) noexcept : h(other.h) { other.h = nullptr; }
    FfiHandle& operator=(FfiHandle&& other) noexcept {
        if (this != &other) {
            rfc_ffi_destroy(h);
            h = other.h;
            other.h = nullptr;
        }
        return *this;
    }

    rfc_ffi* get() const { return h; }
};

// Drop a pending closed 4-point (DIN) cycle so res_raw is only genuinely open
// cycles. Same rule as finalize_res_repeated() in rainflow.c; applied for
// every counting method, including HCM/ASTM, to match rfc()["res_raw"].
static void strip_pending_closed_cycle(std::vector<rfc_ffi_residue_s>& residuum_raw) {
    if (residuum_raw.size() < 4) {
        return;
    }
    const size_t idx = residuum_raw.size() - 4;
    uint32_t A = residuum_raw[idx + 0].cls;
    uint32_t B = residuum_raw[idx + 1].cls;
    uint32_t C = residuum_raw[idx + 2].cls;
    uint32_t D = residuum_raw[idx + 3].cls;

    if (B > C) std::swap(B, C);
    if (A > D) std::swap(A, D);

    if (A <= B && C <= D) {
        residuum_raw.erase(residuum_raw.end() - 3, residuum_raw.end() - 1);
    }
}

static py::array_t<double> values_from_residue_vec(const std::vector<rfc_ffi_residue_s>& raw) {
    py::array_t<double> arr(static_cast<py::ssize_t>(raw.size()));
    auto r = arr.mutable_unchecked<1>();
    for (size_t i = 0; i < raw.size(); i++) {
        r(i) = raw[i].value;
    }
    return arr;
}

// ---------------------------------------------------------------------------
// Error helpers
// ---------------------------------------------------------------------------
static std::runtime_error rfc_error(int32_t code, const char* what) {
    return std::runtime_error(std::string(what) + " (" + rfc_ffi_error_string(code) + ")");
}

static void check_ok(int32_t code, const char* what) {
    if (code != RFC_FFI_OK) {
        throw rfc_error(code, what);
    }
}

static int32_t require_residual_method(int residual_method) {
    if (residual_method < RFC_FFI_RES_NONE || residual_method > RFC_FFI_RES_RP_DIN45667) {
        throw py::value_error("Unknown method for handling residue!");
    }
    return residual_method;
}

static void require_spread_damage(int spread_damage) {
    if (spread_damage < RFC_FFI_SD_NONE || spread_damage > RFC_FFI_SD_TRANSIENT_23c) {
        throw py::value_error("Unknown method for handling damage history!");
    }
}

static void require_lc_method(int lc_method) {
    if (lc_method < RFC_FFI_LC_SLOPES_UP || lc_method > RFC_FFI_LC_FVA) {
        throw py::value_error(
            "Parameter 'lc_method' must be 0, 1, 2 or 3 (LCMethod enumeration, not a flag mask)!");
    }
}

static std::vector<rfc_ffi_residue_s> read_residue(const rfc_ffi* h, const char* what) {
    uint32_t n = 0;
    check_ok(rfc_ffi_residue(h, nullptr, &n), what);
    std::vector<rfc_ffi_residue_s> raw(n);
    if (n > 0) {
        check_ok(rfc_ffi_residue(h, raw.data(), &n), what);
        raw.resize(n);
    }
    return raw;
}

static uint32_t read_class_count(const rfc_ffi* h, const char* what) {
    uint32_t n = 0;
    check_ok(rfc_ffi_class_count(h, &n), what);
    return n;
}

static py::dict dict_from_wl_impaired(const rfc_ffi* h) {
    rfc_ffi_wl_s wl_impaired = {};
    check_ok(rfc_ffi_wl_param_get_impaired(h, &wl_impaired), "Error reading Woehler parameters");
    py::dict wl_out;
    wl_out["sx"] = wl_impaired.sx;
    wl_out["nx"] = wl_impaired.nx;
    wl_out["sd"] = wl_impaired.sd;
    wl_out["nd"] = wl_impaired.nd;
    wl_out["k"]  = wl_impaired.k;
    wl_out["k2"] = wl_impaired.k2;
    wl_out["q"]  = wl_impaired.q;
    wl_out["q2"] = wl_impaired.q2;
    wl_out["omission"] = wl_impaired.omission;
    wl_out["D"] = wl_impaired.D;
    return wl_out;
}

static void set_readonly(py::array_t<double>& arr) {
    arr.attr("setflags")(py::arg("write") = false);
}

static py::array_t<double> frozen_copy(const py::array_t<double>& arr) {
    py::array_t<double> copy = arr.attr("copy")().cast<py::array_t<double>>();
    set_readonly(copy);
    return copy;
}

static py::array_t<double> array_from_res_raw(const rfc_ffi* h) {
    std::vector<rfc_ffi_residue_s> raw = read_residue(h, "Error reading residue");
    strip_pending_closed_cycle(raw);
    return values_from_residue_vec(raw);
}

// ---------------------------------------------------------------------------
// SN-curve ("wl") parameter extraction
// ---------------------------------------------------------------------------
// Mirrors the dict documented in rfc.txt / damage_from_rp.txt:
//   dict(sx=1000, nx=1e7, sd=0, nd=np.inf, k=5, k2=k)
static void check_named_value(double& value, const char* name) {
    const std::string n(name);
    if (n == "sd" || n == "nd" || n == "sx" || n == "nx" || n == "omission") {
        if (value < 0.0) {
            throw py::value_error(std::string("Invalid value for `") + name + "`");
        }
    } else if (n == "k" || n == "k2") {
        value = std::fabs(value);
    }
}

// Parses a `wl` dict into a zero-initialized curve. `extended_def` is set
// true if `sx`, `nx` or `omission` were given explicitly, which selects
// wl_init_any() over wl_init_modified() in the caller.
static void parse_wl_dict(const py::dict& d, rfc_ffi_wl_s& wl, bool& extended_def) {
    extended_def = false;

    for (auto item : d) {
        if (!py::isinstance<py::str>(item.first)) {
            throw py::type_error("Only string keys allowed in `wl`");
        }
        const std::string key = item.first.cast<std::string>();

        double value;
        try {
            value = item.second.cast<double>();
        } catch (const py::cast_error&) {
            throw py::type_error("`" + key + "` must be a numeric type");
        }

        if (key == "sd") { check_named_value(value, "sd"); wl.sd = value; }
        else if (key == "nd") { check_named_value(value, "nd"); wl.nd = value; }
        else if (key == "k") { check_named_value(value, "k"); wl.k = value; }
        else if (key == "k2") { check_named_value(value, "k2"); wl.k2 = value; }
        else if (key == "sx") { check_named_value(value, "sx"); wl.sx = value; extended_def = true; }
        else if (key == "nx") { check_named_value(value, "nx"); wl.nx = value; extended_def = true; }
        else if (key == "omission") { check_named_value(value, "omission"); wl.omission = value; extended_def = true; }
        else {
            throw py::value_error("Wrong key used in wl dict: `" + key + "`");
        }
    }

    wl.k  = std::fabs(wl.k);
    wl.k2 = std::fabs(wl.k2);

    if (wl.k == 0.0) {
        wl.k = 5;
    }
    if (wl.sd == 0.0 && wl.sx == 0.0) {
        wl.sx = 1e3;
    }
    if (wl.nd == 0.0 && wl.nx == 0.0) {
        wl.nx = 1e7;
    }
    if (wl.k2 < wl.k) {
        wl.k2 = wl.k;
    }
    if (wl.sx == 0.0 && wl.sd > 0.0) {
        wl.sx = wl.sd;
        wl.sd = 0.0;
    }
    if (wl.nx == 0.0 && wl.nd > 0.0) {
        wl.nx = wl.nd;
        wl.nd = DBL_MAX;
    }
}

static void init_wl(rfc_ffi* h, bool extended_def, const rfc_ffi_wl_s& wl) {
    if (!extended_def) {
        check_ok(rfc_ffi_wl_init_modified(h, wl.sx, wl.nx, wl.k, wl.k2),
                 "Rainflow initialization error");
    } else {
        check_ok(rfc_ffi_wl_init_any(h, &wl), "Rainflow initialization error");
    }
}

// Flags first, then the level-crossing method. set_lc_method rewrites the
// LC flag bits and stores the conversion mode.
static void apply_flags_and_lc(rfc_ffi* h, int lc_method, bool auto_resize, bool enforce_margin) {
    require_lc_method(lc_method);
    int32_t flags = 0;
    check_ok(rfc_ffi_get_flags(h, &flags), "Rainflow initialization error");
    if (auto_resize) flags |= RFC_FFI_FLAG_AUTORESIZE;
    else             flags &= ~RFC_FFI_FLAG_AUTORESIZE;
    if (enforce_margin) flags |= RFC_FFI_FLAG_ENFORCE_MARGIN;
    else                flags &= ~RFC_FFI_FLAG_ENFORCE_MARGIN;
    check_ok(rfc_ffi_set_flags(h, flags), "Rainflow initialization error");
    check_ok(rfc_ffi_set_lc_method(h, lc_method), "Rainflow initialization error");
}

static void apply_counting_method(rfc_ffi* h, bool use_HCM, bool use_ASTM) {
    if (use_HCM && use_ASTM) {
        throw py::value_error("`use_HCM` and `use_ASTM` are mutually exclusive!");
    }
    if (use_HCM) {
        check_ok(rfc_ffi_set_counting_method(h, RFC_FFI_COUNT_HCM),
                 "Rainflow initialization error");
    }
    if (use_ASTM) {
        check_ok(rfc_ffi_set_counting_method(h, RFC_FFI_COUNT_ASTM),
                 "Rainflow initialization error");
    }
}

static py::array_t<double> array_from_rp(const rfc_ffi* h, const char* what);
static py::array_t<double> array_from_lc(const rfc_ffi* h, const char* what);
static py::array_t<double> array_from_rfm(const rfc_ffi* h, const char* what);
static py::array_t<double> array_from_tp(const rfc_ffi* h);
static py::array_t<double> array_from_residue(const rfc_ffi* h, const char* what);

// ---------------------------------------------------------------------------
// rfc()
// ---------------------------------------------------------------------------
static py::dict rfc(
    py::array_t<double, py::array::c_style | py::array::forcecast> data,
    double class_width,
    int class_count = 100,
    std::optional<double> class_offset = std::nullopt,
    std::optional<double> hysteresis = std::nullopt,
    int residual_method = RFC_FFI_RES_REPEATED,
    int spread_damage = RFC_FFI_SD_TRANSIENT_23c,
    int lc_method = RFC_FFI_LC_SLOPES_ALL,
    bool use_HCM = false,
    bool use_ASTM = false,
    bool enforce_margin = true,
    bool auto_resize = false,
    std::optional<py::dict> wl = std::nullopt
) {
    py::buffer_info buf = data.request();
    if (buf.ndim != 1) {
        throw py::value_error("data must be a 1-D array");
    }
    const auto* ptr = static_cast<const double*>(buf.ptr);
    const size_t len = static_cast<size_t>(buf.shape[0]);

    // class_offset / hysteresis: None means "let the library compute a default",
    // matching the Optional[...] = None semantics in the .pyi (hysteresis
    // defaults to class_width, same as the old `hysteresis < 0` sentinel).
    const double offset = class_offset.value_or(0.0);
    const double hyst    = hysteresis.value_or(class_width);

    rfc_ffi_wl_s wl_param = {};
    wl_param.sx = 1e3;
    wl_param.nx = 1e7;
    wl_param.k = 5;
    wl_param.k2 = 5;
    wl_param.nd = DBL_MAX;
    bool wl_extended_def = false;
    if (wl.has_value()) {
        wl_param = {};
        parse_wl_dict(*wl, wl_param, wl_extended_def);
    }

    FfiHandle rf(rfc_ffi_create());
    if (!rf.get()) {
        throw std::bad_alloc();
    }

    check_ok(rfc_ffi_init(rf.get(), static_cast<uint32_t>(class_count), class_width, offset, hyst,
                          RFC_FFI_FLAGS_DEFAULT),
             "Rainflow initialization error");

    apply_flags_and_lc(rf.get(), lc_method, auto_resize, enforce_margin);
    init_wl(rf.get(), wl_extended_def, wl_param);

    require_spread_damage(spread_damage);
    if (spread_damage > RFC_FFI_SD_NONE) {
        if (rfc_ffi_dh_init(rf.get(), spread_damage, static_cast<uint64_t>(len)) != RFC_FFI_OK) {
            throw std::bad_alloc();
        }
    }

    const int32_t res_method = require_residual_method(residual_method);
    apply_counting_method(rf.get(), use_HCM, use_ASTM);

    std::vector<rfc_ffi_residue_s> residuum_raw;
    {
        // Release the GIL for the actual counting.
        py::gil_scoped_release release;

        check_ok(rfc_ffi_feed(rf.get(), ptr, static_cast<uint64_t>(len)), "Error while counting");
        residuum_raw = read_residue(rf.get(), "Error while counting");
        strip_pending_closed_cycle(residuum_raw);
        check_ok(rfc_ffi_finalize(rf.get(), res_method), "Error while counting");
    }

    // ------------------------------------------------------------------
    // Prepare results
    // ------------------------------------------------------------------
    const char* prep = "Preparing range pair counting";

    py::dict result;

    double damage = 0.0;
    check_ok(rfc_ffi_damage(rf.get(), &damage, nullptr), prep);
    result["damage"] = damage;

    result["rp"] = array_from_rp(rf.get(), prep);
    result["lc"] = array_from_lc(rf.get(), prep);
    result["tp"] = array_from_tp(rf.get());
    result["res_raw"] = values_from_residue_vec(residuum_raw);
    result["res"] = array_from_residue(rf.get(), prep);
    result["rfm"] = array_from_rfm(rf.get(), prep);

    // Damage history, adjacent to the input time series. Missing samples stay 0.
    uint64_t dh_cnt = 0;
    check_ok(rfc_ffi_dh(rf.get(), nullptr, &dh_cnt), prep);
    std::vector<double> dh(static_cast<size_t>(dh_cnt));
    if (dh_cnt > 0) {
        check_ok(rfc_ffi_dh(rf.get(), dh.data(), &dh_cnt), prep);
        dh.resize(static_cast<size_t>(dh_cnt));
    }
    py::array_t<double> dh_arr(static_cast<py::ssize_t>(len));
    {
        auto r = dh_arr.mutable_unchecked<1>();
        for (size_t i = 0; i < len; i++) {
            r(i) = (i < dh.size()) ? dh[i] : 0.0;
        }
    }
    result["dh"] = dh_arr;

    result["wl_miner_consistent"] = dict_from_wl_impaired(rf.get());

    return result;
}

// ---------------------------------------------------------------------------
// damage_from_rp()
// ---------------------------------------------------------------------------
static double damage_from_rp(
    py::array_t<double, py::array::c_style | py::array::forcecast> Sa,
    py::array_t<double, py::array::c_style | py::array::forcecast> counts,
    std::optional<py::dict> wl = std::nullopt,
    int method = 0
) {
    if (method < 0 || method > 3) {
        throw py::value_error("`method` must be in range 0 to 3.");
    }

    rfc_ffi_wl_s wl_param = {};
    if (wl.has_value()) {
        bool extended_def = false;
        parse_wl_dict(*wl, wl_param, extended_def);
    } else {
        wl_param.k  = wl_param.k2 = 5;
        wl_param.sx = 1e3;
        wl_param.nx = 1e7;
        wl_param.sd = 0.0;
        wl_param.nd = DBL_MAX;
    }

    py::buffer_info sa_buf = Sa.request();
    py::buffer_info counts_buf = counts.request();
    if (sa_buf.ndim != 1 || counts_buf.ndim != 1 || sa_buf.shape[0] != counts_buf.shape[0]) {
        throw py::value_error("`Sa` and `counts` must be 1-D arrays of equal length");
    }

    const auto* sa_ptr = static_cast<const double*>(sa_buf.ptr);
    const auto* counts_ptr = static_cast<const double*>(counts_buf.ptr);
    const size_t n = static_cast<size_t>(sa_buf.shape[0]);

    // Sort by amplitude ascending, same as the old ArgSort()-based approach,
    // and carry `counts` along under the same permutation.
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; i++) {
        order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(),
                      [&](size_t a, size_t b) { return sa_ptr[a] < sa_ptr[b]; });

    std::vector<double> vec_sa(n);
    std::vector<double> vec_counts(n);
    for (size_t i = 0; i < n; i++) {
        const size_t src = order[i];
        vec_sa[i] = sa_ptr[src];
        if (counts_ptr[src] < 0) {
            throw py::value_error("Negative values in `counts`.");
        }
        vec_counts[i] = counts_ptr[src];
    }

    FfiHandle rf(rfc_ffi_create());
    if (!rf.get()) {
        throw std::bad_alloc();
    }

    check_ok(rfc_ffi_init(rf.get(), static_cast<uint32_t>(n), 1, 0, 0, RFC_FFI_FLAGS_DEFAULT),
             "Rainflow initialization error");
    check_ok(rfc_ffi_wl_init_any(rf.get(), &wl_param), "Rainflow initialization error");

    double damage = 0.0;
    check_ok(rfc_ffi_damage_from_rp(rf.get(), &damage,
                                    vec_counts.empty() ? nullptr : vec_counts.data(),
                                    vec_sa.empty() ? nullptr : vec_sa.data(),
                                    method),
             "Error while calculation");
    return damage;
}

// ---------------------------------------------------------------------------
// at_transform() — mean-stress correction (Haigh / FKM)
// ---------------------------------------------------------------------------
using RfcArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

struct AtRefCurve {
    std::vector<double> Sa;
    std::vector<double> Sm;

    const double* sa_ptr() const { return Sa.empty() ? nullptr : Sa.data(); }
    const double* sm_ptr() const { return Sm.empty() ? nullptr : Sm.data(); }
    uint32_t count() const { return static_cast<uint32_t>(Sa.size()); }
};

static AtRefCurve parse_at_ref(const std::optional<RfcArray>& Sa_ref,
                               const std::optional<RfcArray>& Sm_ref) {
    AtRefCurve ref;
    if (Sa_ref.has_value() != Sm_ref.has_value()) {
        throw py::value_error("`Sa_ref` and `Sm_ref` must both be given or both omitted");
    }
    if (!Sa_ref.has_value()) {
        return ref;
    }
    const py::buffer_info sa_buf = Sa_ref->request();
    const py::buffer_info sm_buf = Sm_ref->request();
    if (sa_buf.ndim != 1 || sm_buf.ndim != 1 || sa_buf.size != sm_buf.size) {
        throw py::value_error("`Sa_ref` and `Sm_ref` must be 1-D arrays of equal length");
    }
    if (sa_buf.size < 2) {
        throw py::value_error("`Sa_ref` and `Sm_ref` must contain at least 2 points");
    }
    const auto* sa = static_cast<const double*>(sa_buf.ptr);
    const auto* sm = static_cast<const double*>(sm_buf.ptr);
    ref.Sa.assign(sa, sa + sa_buf.size);
    ref.Sm.assign(sm, sm + sm_buf.size);
    return ref;
}

static void at_init_from_args(rfc_ffi* h, double M, double R_rig, double Sm_rig,
                              bool R_pinned, bool symmetric, const AtRefCurve& ref) {
    if (M < 0.0) {
        throw py::value_error("`M` (mean stress sensitivity) must be >= 0");
    }
    if (!ref.Sa.empty() && symmetric) {
        throw py::value_error("`symmetric` is not supported with a custom Haigh reference curve");
    }
    check_ok(rfc_ffi_at_init(h, ref.sa_ptr(), ref.sm_ptr(), ref.count(),
                             M, Sm_rig, R_rig, R_pinned ? 1 : 0, symmetric ? 1 : 0),
             "Amplitude transformation initialization error");
}

static py::array_t<double> at_transform_apply(rfc_ffi* h, const RfcArray& Sa, const RfcArray& Sm) {
    const py::buffer_info sa_buf = Sa.request();
    const py::buffer_info sm_buf = Sm.request();
    if (sa_buf.size != sm_buf.size) {
        throw py::value_error("`Sa` and `Sm` must have the same number of elements");
    }

    py::array_t<double> out(sa_buf.shape);
    const auto* sa_ptr = static_cast<const double*>(sa_buf.ptr);
    const auto* sm_ptr = static_cast<const double*>(sm_buf.ptr);
    double* out_ptr = static_cast<double*>(out.request().ptr);
    int32_t rc = RFC_FFI_OK;
    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < sa_buf.size; i++) {
            rc = rfc_ffi_at_transform(h, sa_ptr[i], sm_ptr[i], &out_ptr[i]);
            if (rc != RFC_FFI_OK) {
                break;
            }
        }
    }
    if (rc != RFC_FFI_OK) {
        throw rfc_error(rc, "Error while amplitude transformation");
    }
    return out;
}

static py::array_t<double> at_transform(
    RfcArray Sa,
    RfcArray Sm,
    double M,
    double R_rig = -1.0,
    double Sm_rig = 0.0,
    bool R_pinned = true,
    std::optional<RfcArray> Sa_ref = std::nullopt,
    std::optional<RfcArray> Sm_ref = std::nullopt,
    bool symmetric = false)
{
    const AtRefCurve ref = parse_at_ref(Sa_ref, Sm_ref);

    FfiHandle rf(rfc_ffi_create());
    if (!rf.get()) {
        throw std::bad_alloc();
    }
    check_ok(rfc_ffi_init(rf.get(), 0, 0.0, 0.0, 0.0, RFC_FFI_FLAGS_DEFAULT),
             "Rainflow initialization error");
    at_init_from_args(rf.get(), M, R_rig, Sm_rig, R_pinned, symmetric, ref);
    return at_transform_apply(rf.get(), Sa, Sm);
}

static py::array_t<double> array_from_rp(const rfc_ffi* h, const char* what) {
    const uint32_t n = read_class_count(h, what);
    std::vector<double> counts(n);
    std::vector<double> sa(n);
    if (n > 0) {
        check_ok(rfc_ffi_rp(h, counts.data(), sa.data()), what);
    }
    py::array_t<double> rp({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(2)});
    auto r = rp.mutable_unchecked<2>();
    for (uint32_t i = 0; i < n; i++) {
        r(i, 0) = sa[i] * 2;  // range = 2 * amplitude
        r(i, 1) = counts[i];
    }
    return rp;
}

static py::array_t<double> array_from_lc(const rfc_ffi* h, const char* what) {
    const uint32_t n = read_class_count(h, what);
    std::vector<double> counts(n);
    std::vector<double> level(n);
    if (n > 0) {
        check_ok(rfc_ffi_lc(h, counts.data(), level.data()), what);
    }
    py::array_t<double> lc({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(2)});
    auto r = lc.mutable_unchecked<2>();
    for (uint32_t i = 0; i < n; i++) {
        r(i, 0) = level[i];
        r(i, 1) = counts[i];
    }
    return lc;
}

static py::array_t<double> array_from_rfm(const rfc_ffi* h, const char* what) {
    const uint32_t n = read_class_count(h, what);
    py::array_t<double> rfm_arr({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(n)});
    if (n == 0) {
        return rfm_arr;
    }
    const size_t cells = static_cast<size_t>(n) * static_cast<size_t>(n);
    std::vector<double> raw(cells);
    check_ok(rfc_ffi_rfm(h, raw.data()), what);
    const double inc = rfc_ffi_full_cycle_increment();
    auto r = rfm_arr.mutable_unchecked<2>();
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            const double counts = raw[static_cast<size_t>(i) * n + j];
            r(i, j) = (inc != 0.0) ? counts / inc : counts;
        }
    }
    return rfm_arr;
}

static py::array_t<double> array_from_tp(const rfc_ffi* h) {
    uint64_t n = 0;
    check_ok(rfc_ffi_tp(h, nullptr, &n), "Error reading turning points");
    std::vector<rfc_ffi_tp_s> tp_buf(static_cast<size_t>(n));
    if (n > 0) {
        check_ok(rfc_ffi_tp(h, tp_buf.data(), &n), "Error reading turning points");
        tp_buf.resize(static_cast<size_t>(n));
    }
    py::array_t<double> tp({static_cast<py::ssize_t>(tp_buf.size()), static_cast<py::ssize_t>(4)});
    auto r = tp.mutable_unchecked<2>();
    for (size_t i = 0; i < tp_buf.size(); i++) {
        r(i, 0) = static_cast<double>(tp_buf[i].pos);
        r(i, 1) = tp_buf[i].value;
        r(i, 2) = tp_buf[i].damage;
        r(i, 3) = static_cast<double>(tp_buf[i].adj_pos);
    }
    return tp;
}

static py::array_t<double> array_from_residue(const rfc_ffi* h, const char* what) {
    const std::vector<rfc_ffi_residue_s> residuum = read_residue(h, what);
    py::array_t<double> arr(static_cast<py::ssize_t>(residuum.size()));
    auto r = arr.mutable_unchecked<1>();
    for (size_t i = 0; i < residuum.size(); i++) {
        r(i) = residuum[i].value;
    }
    return arr;
}

// ---------------------------------------------------------------------------
// Stateful RFC type
// ---------------------------------------------------------------------------
class RFC {
public:
    RFC(double class_width,
        int class_count = 100,
        std::optional<double> class_offset = std::nullopt,
        std::optional<double> hysteresis = std::nullopt,
        bool enforce_margin = true,
        bool auto_resize = false,
        bool use_HCM = false,
        bool use_ASTM = false,
        int spread_damage = RFC_FFI_SD_NONE,
        int lc_method = RFC_FFI_LC_SLOPES_ALL,
        std::optional<py::dict> wl = std::nullopt)
        : rf_(rfc_ffi_create())
    {
        if (!rf_.get()) {
            throw std::bad_alloc();
        }
        configure(class_width, class_count, class_offset, hysteresis,
                  enforce_margin, auto_resize, use_HCM, use_ASTM,
                  spread_damage, lc_method, wl);
    }

    RFC(const RFC&) = delete;
    RFC& operator=(const RFC&) = delete;

    void feed(py::array_t<double, py::array::c_style | py::array::forcecast> data) {
        ensure_feedable();
        py::buffer_info buf = data.request();
        if (buf.ndim != 1) {
            throw py::value_error("data must be a 1-D array");
        }
        const auto* ptr = static_cast<const double*>(buf.ptr);
        const size_t len = static_cast<size_t>(buf.shape[0]);
        int32_t rc = RFC_FFI_OK;
        {
            py::gil_scoped_release release;
            rc = rfc_ffi_feed(rf_.get(), ptr, static_cast<uint64_t>(len));
        }
        check_ok(rc, "Error while counting");
    }

    void finalize(int residual_method = RFC_FFI_RES_REPEATED) {
        ensure_feedable();
        const int32_t res_method = require_residual_method(residual_method);
        res_raw_arr_ = array_from_res_raw(rf_.get());
        set_readonly(res_raw_arr_);
        has_res_raw_ = true;
        int32_t rc = RFC_FFI_OK;
        {
            py::gil_scoped_release release;
            rc = rfc_ffi_finalize(rf_.get(), res_method);
        }
        if (rc != RFC_FFI_OK) {
            has_res_raw_ = false;
            throw rfc_error(rc, "Error while finalizing");
        }
        closed_ = true;
    }

    void close() {
        if (rf_.get() && rfc_ffi_state(rf_.get()) >= RFC_FFI_STATE_INIT) {
            check_ok(rfc_ffi_deinit(rf_.get()), "Error while closing");
        }
        closed_ = true;
    }

    RFC& enter() {
        return *this;
    }

    bool exit(const py::object&, const py::object&, const py::object&) {
        close();
        return false;
    }

    int state() const {
        return rfc_ffi_state(rf_.get());
    }

    int error() const {
        return rfc_ffi_error(rf_.get());
    }

    double damage() const {
        ensure_alive();
        double value = 0.0;
        check_ok(rfc_ffi_damage(rf_.get(), &value, nullptr), "Error reading damage");
        return value;
    }

    py::array_t<double> residue() const {
        ensure_alive();
        return array_from_residue(rf_.get(), "Error reading residue");
    }

    py::array_t<double> rp() const {
        ensure_alive();
        return array_from_rp(rf_.get(), "Error reading range pairs");
    }

    py::array_t<double> lc() const {
        ensure_alive();
        return array_from_lc(rf_.get(), "Error reading level crossings");
    }

    py::array_t<double> rfm() const {
        ensure_alive();
        return array_from_rfm(rf_.get(), "Error reading rainflow matrix");
    }

    py::array_t<double> tp() const {
        ensure_alive();
        return array_from_tp(rf_.get());
    }

    py::array_t<double> res_raw() const {
        ensure_alive();
        if (has_res_raw_) {
            return frozen_copy(res_raw_arr_);
        }
        return array_from_res_raw(rf_.get());
    }

    py::dict wl_miner_consistent() const {
        ensure_alive();
        return dict_from_wl_impaired(rf_.get());
    }

    double damage_as(int residual_method = RFC_FFI_RES_REPEATED) {
        FfiHandle clone = preview_finalized(residual_method);
        double value = 0.0;
        check_ok(rfc_ffi_damage(clone.get(), &value, nullptr), "Error reading preview damage");
        return value;
    }

    py::array_t<double> rp_as(int residual_method = RFC_FFI_RES_REPEATED) {
        return preview_array(residual_method, "Error reading range pairs", array_from_rp);
    }

    py::array_t<double> lc_as(int residual_method = RFC_FFI_RES_REPEATED) {
        return preview_array(residual_method, "Error reading level crossings", array_from_lc);
    }

    py::array_t<double> rfm_as(int residual_method = RFC_FFI_RES_REPEATED) {
        return preview_array(residual_method, "Error reading rainflow matrix", array_from_rfm);
    }

    void at_init(double M,
                 double R_rig = -1.0,
                 double Sm_rig = 0.0,
                 bool R_pinned = true,
                 std::optional<RfcArray> Sa_ref = std::nullopt,
                 std::optional<RfcArray> Sm_ref = std::nullopt,
                 bool symmetric = false)
    {
        ensure_alive();
        if (rfc_ffi_state(rf_.get()) != RFC_FFI_STATE_INIT) {
            throw std::runtime_error(
                "at_init() must be called after construction and before feed()");
        }
        at_ref_ = parse_at_ref(Sa_ref, Sm_ref);
        at_init_from_args(rf_.get(), M, R_rig, Sm_rig, R_pinned, symmetric, at_ref_);
    }

    py::array_t<double> at_transform(RfcArray Sa, RfcArray Sm) {
        ensure_alive();
        return at_transform_apply(rf_.get(), Sa, Sm);
    }

private:
    FfiHandle rf_;
    bool closed_ = false;
    bool has_res_raw_ = false;
    py::array_t<double> res_raw_arr_;
    AtRefCurve at_ref_;

    void ensure_feedable() const {
        if (!rf_.get() || closed_) {
            throw std::runtime_error("RFC object is finalized or closed");
        }
    }

    void ensure_alive() const {
        if (!rf_.get() || rfc_ffi_state(rf_.get()) < RFC_FFI_STATE_INIT) {
            throw std::runtime_error("RFC object is closed");
        }
    }

    FfiHandle preview_finalized(int residual_method) const {
        ensure_feedable();
        const int32_t res_method = require_residual_method(residual_method);
        rfc_ffi* raw = nullptr;
        if (rfc_ffi_clone(rf_.get(), &raw) != RFC_FFI_OK || raw == nullptr) {
            rfc_ffi_destroy(raw);
            throw std::runtime_error("Failed to clone rainflow context");
        }
        FfiHandle clone(raw);
        int32_t rc = RFC_FFI_OK;
        {
            py::gil_scoped_release release;
            rc = rfc_ffi_finalize(clone.get(), res_method);
        }
        check_ok(rc, "Error while preview finalize");
        return clone;
    }

    py::array_t<double> preview_array(int residual_method,
                                      const char* what,
                                      py::array_t<double> (*fn)(const rfc_ffi*, const char*)) {
        FfiHandle clone = preview_finalized(residual_method);
        return fn(clone.get(), what);
    }

    void configure(double class_width,
                   int class_count,
                   std::optional<double> class_offset,
                   std::optional<double> hysteresis,
                   bool enforce_margin,
                   bool auto_resize,
                   bool use_HCM,
                   bool use_ASTM,
                   int spread_damage,
                   int lc_method,
                   const std::optional<py::dict>& wl)
    {
        const double offset = class_offset.value_or(0.0);
        const double hyst = hysteresis.value_or(class_width);

        if (use_HCM && use_ASTM) {
            throw py::value_error("`use_HCM` and `use_ASTM` are mutually exclusive!");
        }
        require_spread_damage(spread_damage);
        // Chunked feed cannot keep a single dh_istream; spread_damage walks
        // that buffer by historical sample position. Use one-shot rfc() instead.
        if (spread_damage > RFC_FFI_SD_NONE) {
            throw py::value_error(
                "Damage history (spread_damage) is not supported on the stateful RFC class; "
                "chunked feed() cannot keep the original sample stream. Use rfc() for "
                "one-shot counting with damage history.");
        }

        rfc_ffi_wl_s wl_param = {};
        wl_param.sx = 1e3;
        wl_param.nx = 1e7;
        wl_param.k = 5;
        wl_param.k2 = 5;
        wl_param.nd = DBL_MAX;
        bool wl_extended_def = false;
        if (wl.has_value()) {
            wl_param = {};
            parse_wl_dict(*wl, wl_param, wl_extended_def);
        }

        check_ok(rfc_ffi_init(rf_.get(), static_cast<uint32_t>(class_count), class_width, offset, hyst,
                              RFC_FFI_FLAGS_DEFAULT),
                 "Rainflow initialization error");
        apply_flags_and_lc(rf_.get(), lc_method, auto_resize, enforce_margin);
        init_wl(rf_.get(), wl_extended_def, wl_param);
        apply_counting_method(rf_.get(), use_HCM, use_ASTM);
    }
};

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------
PYBIND11_MODULE(rfcnt, m) {
    m.doc() = "Python interface for rainflow counting";

    m.def("rfc", &rfc,
          py::arg("data"),
          py::arg("class_width"),
          py::kw_only(),
          py::arg("class_count") = 100,
          py::arg("class_offset") = py::none(),
          py::arg("hysteresis") = py::none(),
          py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED),
          py::arg("spread_damage") = static_cast<int>(RFC_FFI_SD_TRANSIENT_23c),
          py::arg("lc_method") = static_cast<int>(RFC_FFI_LC_SLOPES_ALL),
          py::arg("use_HCM") = false,
          py::arg("use_ASTM") = false,
          py::arg("enforce_margin") = true,
          py::arg("auto_resize") = false,
          py::arg("wl") = py::none(),
          RFC_DOC);

    m.def("damage_from_rp", &damage_from_rp,
          py::arg("Sa"),
          py::arg("counts"),
          py::kw_only(),
          py::arg("wl") = py::none(),
          py::arg("method") = 0,
          DAMAGE_FROM_RP_DOC);

    m.def("at_transform", &at_transform,
          py::arg("Sa"),
          py::arg("Sm"),
          py::kw_only(),
          py::arg("M"),
          py::arg("R_rig") = -1.0,
          py::arg("Sm_rig") = 0.0,
          py::arg("R_pinned") = true,
          py::arg("Sa_ref") = py::none(),
          py::arg("Sm_ref") = py::none(),
          py::arg("symmetric") = false,
          AT_TRANSFORM_DOC);

    py::class_<RFC>(m, "RFC")
        .def(py::init<double, int, std::optional<double>, std::optional<double>,
                      bool, bool, bool, bool, int, int,
                      std::optional<py::dict>>(),
             py::arg("class_width"),
             py::kw_only(),
             py::arg("class_count") = 100,
             py::arg("class_offset") = py::none(),
             py::arg("hysteresis") = py::none(),
             py::arg("enforce_margin") = true,
             py::arg("auto_resize") = false,
             py::arg("use_HCM") = false,
             py::arg("use_ASTM") = false,
             py::arg("spread_damage") = static_cast<int>(RFC_FFI_SD_NONE),
             py::arg("lc_method") = static_cast<int>(RFC_FFI_LC_SLOPES_ALL),
             py::arg("wl") = py::none())
        .def("feed", &RFC::feed, py::arg("data"))
        .def("finalize", &RFC::finalize,
             py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED))
        .def("close", &RFC::close)
        .def("damage_as", &RFC::damage_as,
             py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED))
        .def("rp_as", &RFC::rp_as,
             py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED))
        .def("lc_as", &RFC::lc_as,
             py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED))
        .def("rfm_as", &RFC::rfm_as,
             py::arg("residual_method") = static_cast<int>(RFC_FFI_RES_REPEATED))
        .def("at_init", &RFC::at_init,
             py::arg("M"),
             py::kw_only(),
             py::arg("R_rig") = -1.0,
             py::arg("Sm_rig") = 0.0,
             py::arg("R_pinned") = true,
             py::arg("Sa_ref") = py::none(),
             py::arg("Sm_ref") = py::none(),
             py::arg("symmetric") = false)
        .def("at_transform", &RFC::at_transform,
             py::arg("Sa"),
             py::arg("Sm"))
        .def("__enter__", &RFC::enter, py::return_value_policy::reference_internal)
        .def("__exit__", &RFC::exit)
        .def_property_readonly("state", &RFC::state)
        .def_property_readonly("error", &RFC::error)
        .def_property_readonly("damage", &RFC::damage)
        .def_property_readonly("residue", &RFC::residue)
        .def_property_readonly("rp", &RFC::rp)
        .def_property_readonly("lc", &RFC::lc)
        .def_property_readonly("rfm", &RFC::rfm)
        .def_property_readonly("tp", &RFC::tp)
        .def_property_readonly("res_raw", &RFC::res_raw)
        .def_property_readonly("wl_miner_consistent", &RFC::wl_miner_consistent);
}
