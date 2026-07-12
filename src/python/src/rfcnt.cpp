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
        default:                                     return "Unexpected error";
    }
}

static std::runtime_error rfc_error(const Rainflow& rf, const char* what) {
    return std::runtime_error(std::string(what) + " (" + rfc_err_str(rf.error_get()) + ")");
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
            throw std::runtime_error("Wrong key used in wl dict: `" + key + "`");
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
    int lc_method = 0,
    bool use_HCM = false,
    bool use_ASTM = false,
    bool enforce_margin = true,
    bool auto_resize = false,
    std::optional<py::dict> wl = std::nullopt
) {
    py::buffer_info buf = data.request();
    if (buf.ndim != 1) {
        throw std::runtime_error("data must be a 1-D array");
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

    flags &= ~Rainflow::RFC_FLAGS_COUNT_LC;
    switch (lc_method) {
        case 0: flags |= Rainflow::RFC_FLAGS_COUNT_LC_UP; break;
        case 1: flags |= Rainflow::RFC_FLAGS_COUNT_LC_DN; break;
        case 2: flags |= Rainflow::RFC_FLAGS_COUNT_LC;    break;
        default: throw std::runtime_error("Parameter 'lc_method' must be 0, 1 or 2!");
    }

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

    if (spread_damage < static_cast<int>(Rainflow::RFC_SD_NONE) ||
        spread_damage >= static_cast<int>(Rainflow::RFC_SD_COUNT)) {
        throw std::runtime_error("Unknown method for handling damage history!");
    }
    if (spread_damage > static_cast<int>(Rainflow::RFC_SD_NONE)) {
        if (!rf.dh_init(static_cast<Rainflow::rfc_sd_method_e>(spread_damage), nullptr, len, /*is_static*/ false)) {
            throw std::bad_alloc();
        }
    }

    if (residual_method < static_cast<int>(Rainflow::RFC_RES_NONE) ||
        residual_method >= static_cast<int>(Rainflow::RFC_RES_COUNT)) {
        throw std::runtime_error("Unknown method for handling residue!");
    }

    if (use_HCM && use_ASTM) {
        throw py::value_error("`use_HCM` and `use_ASTM` are mutually exclusive!");
    }
    if (use_HCM)  rf.ctx_get().counting_method = RF::RFC_COUNTING_METHOD_HCM;
    if (use_ASTM) rf.ctx_get().counting_method = RF::RFC_COUNTING_METHOD_ASTM;

    const auto res_method = static_cast<Rainflow::rfc_res_method>(residual_method);
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

        // With regard to finalize_res_repeated() in rainflow.c, remove a
        // pending (already closed) cycle from the raw residuum before it
        // gets reported, so `res_raw` reflects only genuinely open cycles.
        if (residuum_raw.size() >= 4) {
            const size_t idx = residuum_raw.size() - 4;
            unsigned A = residuum_raw[idx + 0].cls;
            unsigned B = residuum_raw[idx + 1].cls;
            unsigned C = residuum_raw[idx + 2].cls;
            unsigned D = residuum_raw[idx + 3].cls;

            if (B > C) std::swap(B, C);
            if (A > D) std::swap(A, D);

            // Check for closed cycles [3]
            if (A <= B && C <= D) {
                residuum_raw.erase(residuum_raw.end() - 3, residuum_raw.end() - 1);
            }
        }

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

    // Turning points: index (1-based), value, pseudo damage, adjacent turning point index.
    auto& tp_storage = rf.tp_storage();
    py::array_t<double> tp({static_cast<py::ssize_t>(tp_storage.size()), static_cast<py::ssize_t>(4)});
    {
        auto r = tp.mutable_unchecked<2>();
        for (size_t i = 0; i < tp_storage.size(); i++) {
            r(i, 0) = static_cast<double>(tp_storage[i].pos);
            r(i, 1) = static_cast<double>(tp_storage[i].value);
            r(i, 2) = static_cast<double>(tp_storage[i].damage);
            r(i, 3) = static_cast<double>(tp_storage[i].adj_pos);
        }
    }
    result["tp"] = tp;

    // Residuum before applying the residual method.
    py::array_t<double> res_raw(static_cast<py::ssize_t>(residuum_raw.size()));
    {
        auto r = res_raw.mutable_unchecked<1>();
        for (size_t i = 0; i < residuum_raw.size(); i++) {
            r(i) = static_cast<double>(residuum_raw[i].value);
        }
    }
    result["res_raw"] = res_raw;

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

    // Miner-consistent (impaired) Woehler curve parameters.
    Rainflow::rfc_wl_param_s wl_impaired;
    if (!rf.wl_param_get_impaired(wl_impaired)) {
        throw rfc_error(rf, "Preparing range pair counting");
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
    result["wl_miner_consistent"] = wl_out;

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
          py::arg("lc_method") = 0,
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
}
