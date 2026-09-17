// rfcnt.cpp — pybind11 port of the rfcnt Python extension.
//
// This wires up the pybind11 interface (argument parsing, defaults, keyword-only
// boundary, GIL handling, return-value construction) to exactly match rfcnt.pyi,
// and ports the algorithm glue (HCM/ASTM variants, spread_damage distribution,
// wl-curve / Miner's-rule handling, and the array-valued outputs: rp, lc, tp,
// res_raw, res, rfm, dh) from the previous raw-CPython-C-API rfcnt.cpp.
//
// Unlike that previous version, we don't talk to the raw C rainflow.h API
// directly: we reuse the same RainflowT<> C++ wrapper (rainflow.hpp) the old
// binding used. It's header-only (no extra translation unit besides
// lib/rainflow.c), and it already owns turning-point storage, Woehler-curve
// bookkeeping and error codes, so re-deriving that against the bare C struct
// would just be re-implementing rainflow.hpp badly.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>   // std::optional <-> None conversion

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <new>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#define RFC_TP_STORAGE std::vector<RF::rfc_value_tuple_s>
#include <rainflow.hpp>

#include "docstrings_generated.hpp"

namespace py = pybind11;

using rfc_residuum_vec = std::vector<Rainflow::rfc_value_tuple_s>;

// Drop a pending closed 4-point (DIN) cycle so res_raw is only genuinely open
// cycles. Same rule as finalize_res_repeated() in rainflow.c; applied for
// every counting method, including HCM/ASTM, to match rfc()["res_raw"].
static void strip_pending_closed_cycle(rfc_residuum_vec& residuum_raw) {
    if (residuum_raw.size() < 4) {
        return;
    }
    const size_t idx = residuum_raw.size() - 4;
    unsigned A = residuum_raw[idx + 0].cls;
    unsigned B = residuum_raw[idx + 1].cls;
    unsigned C = residuum_raw[idx + 2].cls;
    unsigned D = residuum_raw[idx + 3].cls;

    if (B > C) std::swap(B, C);
    if (A > D) std::swap(A, D);

    if (A <= B && C <= D) {
        residuum_raw.erase(residuum_raw.end() - 3, residuum_raw.end() - 1);
    }
}

static py::array_t<double> values_from_residue_vec(const rfc_residuum_vec& raw) {
    py::array_t<double> arr(static_cast<py::ssize_t>(raw.size()));
    auto r = arr.mutable_unchecked<1>();
    for (size_t i = 0; i < raw.size(); i++) {
        r(i) = static_cast<double>(raw[i].value);
    }
    return arr;
}

// ---------------------------------------------------------------------------
// Error helpers
// ---------------------------------------------------------------------------
static const char* rfc_err_str(Rainflow::rfc_error_e nr) {
    switch (nr) {
        case Rainflow::RFC_ERROR_NOERROR:            return "No error";
        case Rainflow::RFC_ERROR_INVARG:             return "Invalid arguments passed";
        case Rainflow::RFC_ERROR_UNSUPPORTED:        return "Unsupported feature";
        case Rainflow::RFC_ERROR_MEMORY:             return "Error on memory allocation";
        case Rainflow::RFC_ERROR_TP:                 return "Error while processing turning points";
        case Rainflow::RFC_ERROR_AT:                 return "Error while amplitude transformation";
        case Rainflow::RFC_ERROR_DH_BAD_STREAM:      return "Input stream must be unique";
        case Rainflow::RFC_ERROR_DH:                 return "Error while damage history calculation/access";
        case Rainflow::RFC_ERROR_LUT:                return "Error while accessing look up tables";
        case Rainflow::RFC_ERROR_DATA_OUT_OF_RANGE:  return "Input data leaves class range";
        case Rainflow::RFC_ERROR_DATA_INCONSISTENT:  return "Processed data is inconsistent";
        default:                                     return "Unexpected error";
    }
}

static std::runtime_error rfc_error(const Rainflow& rf, const char* what) {
    return std::runtime_error(std::string(what) + " (" + rfc_err_str(rf.error_get()) + ")");
}

// lc_method is rfc_lc_count_method (0/1/2/3), not RFC_FLAGS_COUNT_LC_* bits.
// 0..2 = DIN 45667 static slope; 3 = FVA sign-dependent (DIN45667 is a compat alias).
static void apply_lc_method(int lc_method, int& flags, Rainflow& rf) {
    flags &= ~Rainflow::RFC_FLAGS_COUNT_LC;
    switch (lc_method) {
        case Rainflow::RFC_LC_COUNT_METHOD_SLOPES_UP:
            flags |= Rainflow::RFC_FLAGS_COUNT_LC_UP;
            break;
        case Rainflow::RFC_LC_COUNT_METHOD_SLOPES_DOWN:
            flags |= Rainflow::RFC_FLAGS_COUNT_LC_DN;
            break;
        case Rainflow::RFC_LC_COUNT_METHOD_SLOPES_ALL:
            flags |= Rainflow::RFC_FLAGS_COUNT_LC;
            break;
        case Rainflow::RFC_LC_COUNT_METHOD_FVA:
            flags |= Rainflow::RFC_FLAGS_COUNT_LC;
            break;
        default:
            throw py::value_error(
                "Parameter 'lc_method' must be 0, 1, 2 or 3 (LCMethod enumeration, not a flag mask)!");
    }
    rf.ctx_get().lc_count_method =
        static_cast<RF::rfc_lc_count_method_e>(lc_method);
}

static Rainflow::rfc_res_method require_residual_method(int residual_method) {
    if (residual_method < static_cast<int>(Rainflow::RFC_RES_NONE) ||
        residual_method >= static_cast<int>(Rainflow::RFC_RES_COUNT)) {
        throw py::value_error("Unknown method for handling residue!");
    }
    return static_cast<Rainflow::rfc_res_method>(residual_method);
}

static void require_spread_damage(int spread_damage) {
    if (spread_damage < static_cast<int>(Rainflow::RFC_SD_NONE) ||
        spread_damage >= static_cast<int>(Rainflow::RFC_SD_COUNT)) {
        throw py::value_error("Unknown method for handling damage history!");
    }
}

static py::dict dict_from_wl_impaired(Rainflow& rf) {
    Rainflow::rfc_wl_param_s wl_impaired;
    if (!rf.wl_param_get_impaired(wl_impaired)) {
        throw rfc_error(rf, "Error reading Woehler parameters");
    }
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

static py::array_t<double> array_from_res_raw(Rainflow& rf) {
    const Rainflow::rfc_value_tuple_s* residuum = nullptr;
    unsigned residuum_len = 0;
    if (!rf.res_get(&residuum, &residuum_len)) {
        throw rfc_error(rf, "Error reading residue");
    }
    rfc_residuum_vec raw(residuum, residuum + residuum_len);
    strip_pending_closed_cycle(raw);
    return values_from_residue_vec(raw);
}

// Guarantees rf.deinit() runs exactly once, on every exit path (success,
// early throw, or an rf.init() that itself failed) — mirrors the old code's
// unconditional `rf.deinit()` at the bottom of rfc().
struct RainflowDeinitGuard {
    Rainflow* rf;
    ~RainflowDeinitGuard() { rf->deinit(); }
};

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

// Parses a `wl` dict into a zero-initialized rfc_wl_param_s, applying the same
// key validation and cross-field defaulting the old get_dict_wl() used.
// `extended_def` is set true if `sx`, `nx` or `omission` were given explicitly,
// which selects wl_init_any() over wl_init_modified() in the caller.
static void parse_wl_dict(const py::dict& d, Rainflow::rfc_wl_param_s& wl, bool& extended_def) {
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

static py::array_t<double> array_from_tp(Rainflow& rf);

// ---------------------------------------------------------------------------
// rfc()
// ---------------------------------------------------------------------------
static py::dict rfc(
    py::array_t<double, py::array::c_style | py::array::forcecast> data,
    double class_width,
    int class_count = 100,
    std::optional<double> class_offset = std::nullopt,
    std::optional<double> hysteresis = std::nullopt,
    int residual_method = Rainflow::RFC_RES_REPEATED,
    int spread_damage = Rainflow::RFC_SD_TRANSIENT_23c,
    int lc_method = Rainflow::RFC_LC_COUNT_METHOD_SLOPES_ALL,
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

    // Parameters of the SN-curve, if defined.
    double wl_sx = 1e3, wl_nx = 1e7, wl_sd = 0.0, wl_nd = DBL_MAX, wl_k = 5, wl_k2 = 5, wl_omission = 0.0;
    bool wl_extended_def = false;
    if (wl.has_value()) {
        Rainflow::rfc_wl_param_s wl_param = {0};
        parse_wl_dict(*wl, wl_param, wl_extended_def);
        wl_sx = wl_param.sx;
        wl_nx = wl_param.nx;
        wl_sd = wl_param.sd;
        wl_nd = wl_param.nd;
        wl_k  = wl_param.k;
        wl_k2 = wl_param.k2;
        wl_omission = wl_param.omission;
    }

    Rainflow rf;
    RainflowDeinitGuard guard{&rf};

    if (!rf.init(static_cast<unsigned>(class_count), class_width, offset, hyst, Rainflow::RFC_FLAGS_DEFAULT)) {
        throw rfc_error(rf, "Rainflow initialization error");
    }

    int flags = 0;
    rf.flags_get(&flags);

    apply_lc_method(lc_method, flags, rf);

    if (auto_resize) flags |= Rainflow::RFC_FLAGS_AUTORESIZE;
    else             flags &= ~Rainflow::RFC_FLAGS_AUTORESIZE;

    if (enforce_margin) flags |= Rainflow::RFC_FLAGS_ENFORCE_MARGIN;
    else                flags &= ~Rainflow::RFC_FLAGS_ENFORCE_MARGIN;

    rf.flags_set(flags, /* debugging */ false, /* overwrite */ true);

    if (!wl_extended_def) {
        if (!rf.wl_init_modified(wl_sx, wl_nx, wl_k, wl_k2)) {
            throw rfc_error(rf, "Rainflow initialization error");
        }
    } else {
        Rainflow::rfc_wl_param_s wl_param = {0};
        wl_param.sd = wl_sd;
        wl_param.nd = wl_nd;
        wl_param.sx = wl_sx;
        wl_param.nx = wl_nx;
        wl_param.k  = wl_k;
        wl_param.k2 = wl_k2;
        wl_param.omission = wl_omission;
        if (!rf.wl_init_any(&wl_param)) {
            throw rfc_error(rf, "Rainflow initialization error");
        }
    }

    require_spread_damage(spread_damage);
    if (spread_damage > static_cast<int>(Rainflow::RFC_SD_NONE)) {
        if (!rf.dh_init(static_cast<Rainflow::rfc_sd_method_e>(spread_damage), nullptr, len, /*is_static*/ false)) {
            throw std::bad_alloc();
        }
    }

    const auto res_method = require_residual_method(residual_method);

    if (use_HCM && use_ASTM) {
        throw py::value_error("`use_HCM` and `use_ASTM` are mutually exclusive!");
    }
    if (use_HCM)  rf.ctx_get().counting_method = RF::RFC_COUNTING_METHOD_HCM;
    if (use_ASTM) rf.ctx_get().counting_method = RF::RFC_COUNTING_METHOD_ASTM;

    rfc_residuum_vec residuum_raw;
    {
        // Release the GIL for the actual counting — this is the concrete
        // win over the raw C-API version, which held the GIL throughout.
        py::gil_scoped_release release;

        if (!rf.feed(ptr, len)) {
            throw rfc_error(rf, "Error while counting");
        }

        const Rainflow::rfc_value_tuple_s* residuum;
        unsigned residuum_len;
        if (!rf.res_get(&residuum, &residuum_len)) {
            throw rfc_error(rf, "Error while counting");
        }
        residuum_raw.assign(residuum, residuum + residuum_len);
        strip_pending_closed_cycle(residuum_raw);

        if (!rf.finalize(res_method)) {
            throw rfc_error(rf, "Error while counting");
        }
    }

    // ------------------------------------------------------------------
    // Prepare results
    // ------------------------------------------------------------------
    unsigned class_count_actual;
    if (!rf.class_count(&class_count_actual)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }

    py::dict result;

    double damage = 0.0;
    if (!rf.damage(&damage)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    result["damage"] = damage;

    // Range pair counts: first column is range (= 2 * amplitude), ascending.
    Rainflow::rfc_counts_v ct;
    Rainflow::rfc_value_v sa;
    if (!rf.rp_get(ct, sa)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    py::array_t<double> rp({static_cast<py::ssize_t>(class_count_actual), static_cast<py::ssize_t>(2)});
    {
        auto r = rp.mutable_unchecked<2>();
        for (unsigned i = 0; i < class_count_actual; i++) {
            r(i, 0) = static_cast<double>(sa[i]) * 2;  // range = 2 * amplitude
            r(i, 1) = static_cast<double>(ct[i]);
        }
    }
    result["rp"] = rp;

    // Level crossings: first column is class upper limit, ascending.
    if (!rf.lc_get(ct, sa)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    py::array_t<double> lc({static_cast<py::ssize_t>(class_count_actual), static_cast<py::ssize_t>(2)});
    {
        auto r = lc.mutable_unchecked<2>();
        for (unsigned i = 0; i < class_count_actual; i++) {
            r(i, 0) = static_cast<double>(sa[i]);
            r(i, 1) = static_cast<double>(ct[i]);
        }
    }
    result["lc"] = lc;

    result["tp"] = array_from_tp(rf);

    result["res_raw"] = values_from_residue_vec(residuum_raw);

    // Residuum after applying the residual method.
    const Rainflow::rfc_value_tuple_s* p_residue;
    unsigned residue_cnt;
    if (!rf.res_get(&p_residue, &residue_cnt)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    py::array_t<double> res(static_cast<py::ssize_t>(residue_cnt));
    {
        auto r = res.mutable_unchecked<1>();
        for (unsigned i = 0; i < residue_cnt; i++) {
            r(i) = static_cast<double>(p_residue[i].value);
        }
    }
    result["res"] = res;

    // Rainflow matrix, class_count x class_count, counts as full-cycle units.
    Rainflow::rfc_rfm_item_v rfm;
    if (!rf.rfm_get(rfm)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    py::array_t<double> rfm_arr({static_cast<py::ssize_t>(class_count_actual), static_cast<py::ssize_t>(class_count_actual)});
    {
        auto r = rfm_arr.mutable_unchecked<2>();
        for (unsigned i = 0; i < class_count_actual; i++) {
            for (unsigned j = 0; j < class_count_actual; j++) {
                r(i, j) = 0.0;
            }
        }
        for (const auto& item : rfm) {
            r(item.from, item.to) += static_cast<double>(item.counts) / RFC_FULL_CYCLE_INCREMENT;
        }
    }
    result["rfm"] = rfm_arr;

    // Damage history, adjacent to the input time series.
    const double* dh;
    size_t dh_cnt;
    if (!rf.dh_get(&dh, &dh_cnt)) {
        throw rfc_error(rf, "Preparing range pair counting");
    }
    py::array_t<double> dh_arr(static_cast<py::ssize_t>(len));
    {
        auto r = dh_arr.mutable_unchecked<1>();
        for (size_t i = 0; i < len; i++) {
            r(i) = (i < dh_cnt) ? dh[i] : 0.0;
        }
    }
    result["dh"] = dh_arr;

    result["wl_miner_consistent"] = dict_from_wl_impaired(rf);

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

    Rainflow::rfc_wl_param_s wl_param = {0};
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
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                      [&](size_t a, size_t b) { return sa_ptr[a] < sa_ptr[b]; });

    Rainflow::rfc_value_v vec_sa(n);
    Rainflow::rfc_counts_v vec_counts(n);
    for (size_t i = 0; i < n; i++) {
        const size_t src = order[i];
        vec_sa[i] = sa_ptr[src];
        if (counts_ptr[src] < 0) {
            throw py::value_error("Negative values in `counts`.");
        }
        vec_counts[i] = static_cast<Rainflow::rfc_counts_t>(counts_ptr[src]);
    }

    Rainflow rf;
    RainflowDeinitGuard guard{&rf};

    if (!rf.init(static_cast<unsigned>(n), 1, 0, 0)) {
        throw rfc_error(rf, "Rainflow initialization error");
    }
    if (!rf.wl_init_any(&wl_param)) {
        throw rfc_error(rf, "Rainflow initialization error");
    }

    double damage;
    const auto rp_calc_type = static_cast<Rainflow::rfc_rp_damage_method>(method);
    if (!rf.damage_from_rp(damage, vec_counts, vec_sa, rp_calc_type)) {
        throw rfc_error(rf, "Error while calculation");
    }

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
    unsigned count() const { return static_cast<unsigned>(Sa.size()); }
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

static void at_init_from_args(Rainflow& rf, double M, double R_rig, double Sm_rig,
                              bool R_pinned, bool symmetric, const AtRefCurve& ref) {
    if (M < 0.0) {
        throw py::value_error("`M` (mean stress sensitivity) must be >= 0");
    }
    if (!ref.Sa.empty() && symmetric) {
        throw py::value_error("`symmetric` is not supported with a custom Haigh reference curve");
    }
    if (!rf.at_init(ref.sa_ptr(), ref.sm_ptr(), ref.count(),
                    M, Sm_rig, R_rig, R_pinned, symmetric)) {
        throw rfc_error(rf, "Amplitude transformation initialization error");
    }
}

static py::array_t<double> at_transform_apply(Rainflow& rf, const RfcArray& Sa, const RfcArray& Sm) {
    const py::buffer_info sa_buf = Sa.request();
    const py::buffer_info sm_buf = Sm.request();
    if (sa_buf.size != sm_buf.size) {
        throw py::value_error("`Sa` and `Sm` must have the same number of elements");
    }

    py::array_t<double> out(sa_buf.shape);
    const auto* sa_ptr = static_cast<const double*>(sa_buf.ptr);
    const auto* sm_ptr = static_cast<const double*>(sm_buf.ptr);
    double* out_ptr = static_cast<double*>(out.request().ptr);
    for (py::ssize_t i = 0; i < sa_buf.size; i++) {
        double Sa_t = 0.0;
        if (!rf.at_transform(sa_ptr[i], sm_ptr[i], Sa_t)) {
            throw rfc_error(rf, "Error while amplitude transformation");
        }
        out_ptr[i] = Sa_t;
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

    Rainflow rf;
    RainflowDeinitGuard guard{&rf};
    if (!rf.init(0, 0.0, 0.0, 0.0, Rainflow::RFC_FLAGS_DEFAULT)) {
        throw rfc_error(rf, "Rainflow initialization error");
    }
    at_init_from_args(rf, M, R_rig, Sm_rig, R_pinned, symmetric, ref);
    return at_transform_apply(rf, Sa, Sm);
}

// ---------------------------------------------------------------------------
// Heap Rainflow — stack instances have been unsafe with MinGW + CPython.
// ---------------------------------------------------------------------------
static Rainflow* rfc_rainflow_new() {
    void* mem = std::malloc(sizeof(Rainflow));
    if (!mem) {
        return nullptr;
    }
    return new (mem) Rainflow();
}

static void rfc_rainflow_delete(Rainflow* rf) {
    if (!rf) {
        return;
    }
    rf->deinit();
    rf->ctx_get().internal.obj = nullptr;
    rf->~Rainflow();
    std::free(rf);
}

// Deep-copy counting state so a preview finalize can run on the copy.
static bool rfc_rainflow_clone(const Rainflow* src, Rainflow* dst) {
    unsigned class_count = 0;
    Rainflow::rfc_value_t class_width = 0;
    Rainflow::rfc_value_t class_offset = 0;
    Rainflow::rfc_value_t hysteresis = 0;
    int flags = 0;
    Rainflow::rfc_wl_param_s wl = {};

    if (!src || !dst) {
        return false;
    }

    if (!src->class_count(&class_count) ||
        !src->class_width(&class_width) ||
        !src->class_offset(&class_offset) ||
        !src->hysteresis(&hysteresis) ||
        !src->flags_get(&flags)) {
        return false;
    }

    flags &= ~Rainflow::RFC_FLAGS_COUNT_DH;

    if (!dst->init(class_count, class_width, class_offset, hysteresis,
                   static_cast<Rainflow::rfc_flags_e>(flags))) {
        return false;
    }

    if (!src->wl_param_get(wl) || !dst->wl_init_any(&wl)) {
        return false;
    }

    const Rainflow::rfc_ctx_s& s = src->ctx_get();
    Rainflow::rfc_ctx_s& d = dst->ctx_get();

    d.counting_method      = s.counting_method;
    d.lc_count_method      = s.lc_count_method;
    d.series_start         = s.series_start;
    d.series_end           = s.series_end;
    d.series_bounds_valid  = s.series_bounds_valid;
    d.full_inc        = s.full_inc;
    d.half_inc        = s.half_inc;
    d.curr_inc        = s.curr_inc;
    d.damage          = s.damage;
    d.damage_residue  = s.damage_residue;

    if (s.rfm && d.rfm) {
        std::memcpy(d.rfm, s.rfm,
                    static_cast<size_t>(class_count) * static_cast<size_t>(class_count) * sizeof(*d.rfm));
    }
    if (s.rp && d.rp) {
        std::memcpy(d.rp, s.rp, static_cast<size_t>(class_count) * sizeof(*d.rp));
    }
    if (s.lc && d.lc) {
        std::memcpy(d.lc, s.lc, static_cast<size_t>(class_count) * sizeof(*d.lc));
    }

    const size_t res_n = s.residue_cnt + (s.state == RF::RFC_STATE_BUSY_INTERIM ? 1u : 0u);
    if (res_n > d.residue_cap) {
        return false;
    }
    if (res_n && s.residue && d.residue) {
        std::memcpy(d.residue, s.residue, res_n * sizeof(*d.residue));
    }
    d.residue_cnt = s.residue_cnt;

    d.internal.slope      = s.internal.slope;
    d.internal.extrema[0] = s.internal.extrema[0];
    d.internal.extrema[1] = s.internal.extrema[1];
#if RFC_GLOBAL_EXTREMA
    d.internal.extrema_changed = s.internal.extrema_changed;
#endif
    d.internal.pos        = s.internal.pos;
    d.internal.pos_offset = s.internal.pos_offset;
#if RFC_TP_SUPPORT
    d.internal.margin[0]    = s.internal.margin[0];
    d.internal.margin[1]    = s.internal.margin[1];
    d.internal.margin_stage = s.internal.margin_stage;
#endif
#if !RFC_MINIMAL
    d.internal.wl = s.internal.wl;
#endif

#if RFC_AT_SUPPORT
    d.internal.at_haigh = s.internal.at_haigh;
    d.at = s.at;
    if (s.at.Sa == s.internal.at_haigh.Sa) {
        d.at.Sa = d.internal.at_haigh.Sa;
        d.at.Sm = d.internal.at_haigh.Sm;
    }
#if RFC_DAMAGE_FAST
    // LUT was built by dst->init() without AT; force the slow path that
    // calls RFC_at_transform until a later at_init rebuilds the tables.
    d.damage_lut_inapt = 1;
#endif
#endif

#if RFC_HCM_SUPPORT
    d.internal.hcm.IR = s.internal.hcm.IR;
    d.internal.hcm.IZ = s.internal.hcm.IZ;
    if (s.internal.hcm.stack && s.internal.hcm.stack_cap) {
        if (!d.internal.hcm.stack ||
            d.internal.hcm.stack_cap < s.internal.hcm.stack_cap) {
            return false;
        }
        std::memcpy(d.internal.hcm.stack, s.internal.hcm.stack,
                    s.internal.hcm.stack_cap * sizeof(*d.internal.hcm.stack));
    }
#endif

#if RFC_TP_SUPPORT
    dst->tp_storage()    = src->tp_storage();
    d.tp_cnt             = s.tp_cnt;
    d.tp_cap             = dst->tp_storage().capacity();
    d.tp_locked          = s.tp_locked;
    d.tp_prune_size      = s.tp_prune_size;
    d.tp_prune_threshold = s.tp_prune_threshold;
#endif

    (void)dst->flags_unset(Rainflow::RFC_FLAGS_COUNT_DH, /*debugging*/ false);

    d.state = s.state;
    d.error = s.error;
    return true;
}

static py::array_t<double> array_from_rp(Rainflow& rf) {
    unsigned class_count = 0;
    Rainflow::rfc_counts_v ct;
    Rainflow::rfc_value_v sa;
    if (!rf.class_count(&class_count) || !rf.rp_get(ct, sa)) {
        throw rfc_error(rf, "Error reading range pairs");
    }
    py::array_t<double> rp({static_cast<py::ssize_t>(class_count), static_cast<py::ssize_t>(2)});
    auto r = rp.mutable_unchecked<2>();
    for (unsigned i = 0; i < class_count; i++) {
        r(i, 0) = static_cast<double>(sa[i]) * 2;
        r(i, 1) = static_cast<double>(ct[i]);
    }
    return rp;
}

static py::array_t<double> array_from_lc(Rainflow& rf) {
    unsigned class_count = 0;
    Rainflow::rfc_counts_v ct;
    Rainflow::rfc_value_v sa;
    if (!rf.class_count(&class_count) || !rf.lc_get(ct, sa)) {
        throw rfc_error(rf, "Error reading level crossings");
    }
    py::array_t<double> lc({static_cast<py::ssize_t>(class_count), static_cast<py::ssize_t>(2)});
    auto r = lc.mutable_unchecked<2>();
    for (unsigned i = 0; i < class_count; i++) {
        r(i, 0) = static_cast<double>(sa[i]);
        r(i, 1) = static_cast<double>(ct[i]);
    }
    return lc;
}

static py::array_t<double> array_from_rfm(Rainflow& rf) {
    unsigned class_count = 0;
    Rainflow::rfc_rfm_item_v rfm;
    if (!rf.class_count(&class_count) || !rf.rfm_get(rfm)) {
        throw rfc_error(rf, "Error reading rainflow matrix");
    }
    py::array_t<double> rfm_arr({static_cast<py::ssize_t>(class_count),
                                 static_cast<py::ssize_t>(class_count)});
    auto r = rfm_arr.mutable_unchecked<2>();
    for (unsigned i = 0; i < class_count; i++) {
        for (unsigned j = 0; j < class_count; j++) {
            r(i, j) = 0.0;
        }
    }
    for (const auto& item : rfm) {
        r(item.from, item.to) += static_cast<double>(item.counts) / RFC_FULL_CYCLE_INCREMENT;
    }
    return rfm_arr;
}

static py::array_t<double> array_from_tp(Rainflow& rf) {
    const auto& tp_storage = rf.tp_storage();
    size_t n = rf.ctx_get().tp_cnt;
    if (n > tp_storage.size()) {
        n = tp_storage.size();
    }
    py::array_t<double> tp({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(4)});
    auto r = tp.mutable_unchecked<2>();
    for (size_t i = 0; i < n; i++) {
        r(i, 0) = static_cast<double>(tp_storage[i].pos);
        r(i, 1) = static_cast<double>(tp_storage[i].value);
        r(i, 2) = static_cast<double>(tp_storage[i].damage);
        r(i, 3) = static_cast<double>(tp_storage[i].adj_pos);
    }
    return tp;
}

static py::array_t<double> array_from_residue(Rainflow& rf) {
    const Rainflow::rfc_value_tuple_s* residuum = nullptr;
    unsigned residuum_len = 0;
    if (!rf.res_get(&residuum, &residuum_len)) {
        throw rfc_error(rf, "Error reading residue");
    }
    py::array_t<double> arr(static_cast<py::ssize_t>(residuum_len));
    auto r = arr.mutable_unchecked<1>();
    for (unsigned i = 0; i < residuum_len; i++) {
        r(i) = static_cast<double>(residuum[i].value);
    }
    return arr;
}

// ---------------------------------------------------------------------------
// Stateful RFC type (owns a heap RainflowT)
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
        int spread_damage = Rainflow::RFC_SD_NONE,
        int lc_method = Rainflow::RFC_LC_COUNT_METHOD_SLOPES_ALL,
        std::optional<py::dict> wl = std::nullopt)
        : rf_(rfc_rainflow_new())
    {
        if (!rf_) {
            throw std::bad_alloc();
        }
        try {
            configure(class_width, class_count, class_offset, hysteresis,
                      enforce_margin, auto_resize, use_HCM, use_ASTM,
                      spread_damage, lc_method, wl);
        } catch (...) {
            rfc_rainflow_delete(rf_);
            rf_ = nullptr;
            throw;
        }
    }

    RFC(const RFC&) = delete;
    RFC& operator=(const RFC&) = delete;

    ~RFC() {
        rfc_rainflow_delete(rf_);
        rf_ = nullptr;
    }

    void feed(py::array_t<double, py::array::c_style | py::array::forcecast> data) {
        ensure_feedable();
        py::buffer_info buf = data.request();
        if (buf.ndim != 1) {
            throw py::value_error("data must be a 1-D array");
        }
        const auto* ptr = static_cast<const double*>(buf.ptr);
        const size_t len = static_cast<size_t>(buf.shape[0]);
        bool ok = false;
        {
            py::gil_scoped_release release;
            ok = rf_->feed(ptr, len);
        }
        if (!ok) {
            throw rfc_error(*rf_, "Error while counting");
        }
    }

    void finalize(int residual_method = Rainflow::RFC_RES_REPEATED) {
        ensure_feedable();
        const auto res_method = require_residual_method(residual_method);
        res_raw_arr_ = array_from_res_raw(*rf_);
        set_readonly(res_raw_arr_);
        has_res_raw_ = true;
        bool ok = false;
        {
            py::gil_scoped_release release;
            ok = rf_->finalize(res_method);
        }
        if (!ok) {
            has_res_raw_ = false;
            throw rfc_error(*rf_, "Error while finalizing");
        }
        closed_ = true;
    }

    void close() {
        if (rf_ && rf_->state_get() >= Rainflow::RFC_STATE_INIT) {
            if (!rf_->deinit()) {
                throw rfc_error(*rf_, "Error while closing");
            }
            rf_->ctx_get().internal.obj = nullptr;
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
        return static_cast<int>(rf_->state_get());
    }

    int error() const {
        return static_cast<int>(rf_->error_get());
    }

    double damage() const {
        ensure_alive();
        Rainflow::rfc_value_t value = 0;
        if (!rf_->damage(&value, nullptr)) {
            throw rfc_error(*rf_, "Error reading damage");
        }
        return static_cast<double>(value);
    }

    py::array_t<double> residue() const {
        ensure_alive();
        return array_from_residue(*rf_);
    }

    py::array_t<double> rp() const {
        ensure_alive();
        return array_from_rp(*rf_);
    }

    py::array_t<double> lc() const {
        ensure_alive();
        return array_from_lc(*rf_);
    }

    py::array_t<double> rfm() const {
        ensure_alive();
        return array_from_rfm(*rf_);
    }

    py::array_t<double> tp() const {
        ensure_alive();
        return array_from_tp(*rf_);
    }

    py::array_t<double> res_raw() const {
        ensure_alive();
        if (has_res_raw_) {
            return frozen_copy(res_raw_arr_);
        }
        return array_from_res_raw(*rf_);
    }

    py::dict wl_miner_consistent() const {
        ensure_alive();
        return dict_from_wl_impaired(*rf_);
    }

    double damage_as(int residual_method = Rainflow::RFC_RES_REPEATED) {
        Rainflow* clone = preview_finalized(residual_method);
        Rainflow::rfc_value_t value = 0;
        const bool ok = clone->damage(&value, nullptr);
        if (!ok) {
            const auto err = rfc_error(*clone, "Error reading preview damage");
            rfc_rainflow_delete(clone);
            throw err;
        }
        rfc_rainflow_delete(clone);
        return static_cast<double>(value);
    }

    py::array_t<double> rp_as(int residual_method = Rainflow::RFC_RES_REPEATED) {
        return preview_array(residual_method, array_from_rp);
    }

    py::array_t<double> lc_as(int residual_method = Rainflow::RFC_RES_REPEATED) {
        return preview_array(residual_method, array_from_lc);
    }

    py::array_t<double> rfm_as(int residual_method = Rainflow::RFC_RES_REPEATED) {
        return preview_array(residual_method, array_from_rfm);
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
        if (rf_->state_get() != Rainflow::RFC_STATE_INIT) {
            throw std::runtime_error(
                "at_init() must be called after construction and before feed()");
        }
        at_ref_ = parse_at_ref(Sa_ref, Sm_ref);
        at_init_from_args(*rf_, M, R_rig, Sm_rig, R_pinned, symmetric, at_ref_);
    }

    py::array_t<double> at_transform(RfcArray Sa, RfcArray Sm) {
        ensure_alive();
        return at_transform_apply(*rf_, Sa, Sm);
    }

private:
    Rainflow* rf_ = nullptr;
    bool closed_ = false;
    bool has_res_raw_ = false;
    py::array_t<double> res_raw_arr_;
    AtRefCurve at_ref_;

    void ensure_feedable() const {
        if (!rf_ || closed_) {
            throw std::runtime_error("RFC object is finalized or closed");
        }
    }

    void ensure_alive() const {
        if (!rf_ || rf_->state_get() < Rainflow::RFC_STATE_INIT) {
            throw std::runtime_error("RFC object is closed");
        }
    }

    Rainflow* preview_finalized(int residual_method) const {
        ensure_feedable();
        const auto res_method = require_residual_method(residual_method);
        Rainflow* clone = rfc_rainflow_new();
        if (!clone) {
            throw std::bad_alloc();
        }
        if (!rfc_rainflow_clone(rf_, clone)) {
            rfc_rainflow_delete(clone);
            throw std::runtime_error("Failed to clone rainflow context");
        }
        bool ok = false;
        {
            py::gil_scoped_release release;
            ok = clone->finalize(res_method);
        }
        if (!ok) {
            const auto err = rfc_error(*clone, "Error while preview finalize");
            rfc_rainflow_delete(clone);
            throw err;
        }
        return clone;
    }

    py::array_t<double> preview_array(int residual_method,
                                      py::array_t<double> (*fn)(Rainflow&)) {
        Rainflow* clone = preview_finalized(residual_method);
        try {
            py::array_t<double> arr = fn(*clone);
            rfc_rainflow_delete(clone);
            return arr;
        } catch (...) {
            rfc_rainflow_delete(clone);
            throw;
        }
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
        if (spread_damage > static_cast<int>(Rainflow::RFC_SD_NONE)) {
            throw py::value_error(
                "Damage history (spread_damage) is not supported on the stateful RFC class; "
                "chunked feed() cannot keep the original sample stream. Use rfc() for "
                "one-shot counting with damage history.");
        }

        double wl_sx = 1e3, wl_nx = 1e7, wl_sd = 0.0, wl_nd = DBL_MAX;
        double wl_k = 5, wl_k2 = 5, wl_omission = 0.0;
        bool wl_extended_def = false;
        if (wl.has_value()) {
            Rainflow::rfc_wl_param_s wl_param = {0};
            parse_wl_dict(*wl, wl_param, wl_extended_def);
            wl_sx = wl_param.sx;
            wl_nx = wl_param.nx;
            wl_sd = wl_param.sd;
            wl_nd = wl_param.nd;
            wl_k = wl_param.k;
            wl_k2 = wl_param.k2;
            wl_omission = wl_param.omission;
        }

        if (!rf_->init(static_cast<unsigned>(class_count), class_width, offset, hyst,
                       Rainflow::RFC_FLAGS_DEFAULT)) {
            throw rfc_error(*rf_, "Rainflow initialization error");
        }

        int flags = 0;
        rf_->flags_get(&flags);
        apply_lc_method(lc_method, flags, *rf_);
        if (auto_resize) flags |= Rainflow::RFC_FLAGS_AUTORESIZE;
        else             flags &= ~Rainflow::RFC_FLAGS_AUTORESIZE;
        if (enforce_margin) flags |= Rainflow::RFC_FLAGS_ENFORCE_MARGIN;
        else                flags &= ~Rainflow::RFC_FLAGS_ENFORCE_MARGIN;
        rf_->flags_set(flags, /*debugging*/ false, /*overwrite*/ true);

        if (!wl_extended_def) {
            if (!rf_->wl_init_modified(wl_sx, wl_nx, wl_k, wl_k2)) {
                throw rfc_error(*rf_, "Rainflow initialization error");
            }
        } else {
            Rainflow::rfc_wl_param_s wl_param = {0};
            wl_param.sd = wl_sd;
            wl_param.nd = wl_nd;
            wl_param.sx = wl_sx;
            wl_param.nx = wl_nx;
            wl_param.k = wl_k;
            wl_param.k2 = wl_k2;
            wl_param.omission = wl_omission;
            if (!rf_->wl_init_any(&wl_param)) {
                throw rfc_error(*rf_, "Rainflow initialization error");
            }
        }

        if (use_HCM) {
            rf_->ctx_get().counting_method = RF::RFC_COUNTING_METHOD_HCM;
        }
        if (use_ASTM) {
            rf_->ctx_get().counting_method = RF::RFC_COUNTING_METHOD_ASTM;
        }
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
          py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED),
          py::arg("spread_damage") = static_cast<int>(Rainflow::RFC_SD_TRANSIENT_23c),
          py::arg("lc_method") = static_cast<int>(Rainflow::RFC_LC_COUNT_METHOD_SLOPES_ALL),
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
             py::arg("spread_damage") = static_cast<int>(Rainflow::RFC_SD_NONE),
             py::arg("lc_method") = static_cast<int>(Rainflow::RFC_LC_COUNT_METHOD_SLOPES_ALL),
             py::arg("wl") = py::none())
        .def("feed", &RFC::feed, py::arg("data"))
        .def("finalize", &RFC::finalize,
             py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED))
        .def("close", &RFC::close)
        .def("damage_as", &RFC::damage_as,
             py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED))
        .def("rp_as", &RFC::rp_as,
             py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED))
        .def("lc_as", &RFC::lc_as,
             py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED))
        .def("rfm_as", &RFC::rfm_as,
             py::arg("residual_method") = static_cast<int>(Rainflow::RFC_RES_REPEATED))
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
