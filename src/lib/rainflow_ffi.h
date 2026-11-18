/*
 * Stable FFI for the rainflow core.
 *
 * This header does not include rainflow.h and has no feature #ifdefs.
 * The context layout, packing, and compile-time flags stay inside the
 * library. Callers see an opaque handle, fixed-width structs, and an
 * int32_t error code. A feature compiled out of the core returns
 * RFC_FFI_ERR_UNSUPPORTED.
 *
 * Integer values match the full-build enums in rainflow.h.
 */
#ifndef RAINFLOW_FFI_H
#define RAINFLOW_FFI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rfc_ffi rfc_ffi;

/* Error codes. Explicit values, including those a minimal build omits. */
#define RFC_FFI_OK                      0
#define RFC_FFI_ERR_UNEXP              -1
#define RFC_FFI_ERR_INVARG              1
#define RFC_FFI_ERR_UNSUPPORTED         2
#define RFC_FFI_ERR_MEMORY              3
#define RFC_FFI_ERR_TP                  4
#define RFC_FFI_ERR_AT                  5
#define RFC_FFI_ERR_DH_BAD_STREAM       6
#define RFC_FFI_ERR_DH                  7
#define RFC_FFI_ERR_LUT                 8
#define RFC_FFI_ERR_DATA_OUT_OF_RANGE   9
#define RFC_FFI_ERR_DATA_INCONSISTENT  10

/* Residual methods (rfc_res_method). */
#define RFC_FFI_RES_NONE                0
#define RFC_FFI_RES_IGNORE              1
#define RFC_FFI_RES_NO_FINALIZE         2
#define RFC_FFI_RES_DISCARD             3
#define RFC_FFI_RES_HALFCYCLES          4
#define RFC_FFI_RES_FULLCYCLES          5
#define RFC_FFI_RES_CLORMANN_SEEGER     6
#define RFC_FFI_RES_REPEATED            7
#define RFC_FFI_RES_RP_DIN45667         8

/* Spread-damage methods (rfc_sd_method). */
#define RFC_FFI_SD_NONE                -1
#define RFC_FFI_SD_HALF_23              0
#define RFC_FFI_SD_RAMP_AMPLITUDE_23    1
#define RFC_FFI_SD_RAMP_DAMAGE_23       2
#define RFC_FFI_SD_RAMP_AMPLITUDE_24    3
#define RFC_FFI_SD_RAMP_DAMAGE_24       4
#define RFC_FFI_SD_FULL_P2              5
#define RFC_FFI_SD_FULL_P3              6
#define RFC_FFI_SD_TRANSIENT_23         7
#define RFC_FFI_SD_TRANSIENT_23c        8

/* Level-crossing methods (rfc_lc_count_method). 3 is FVA; DIN45667 aliases it. */
#define RFC_FFI_LC_SLOPES_UP            0
#define RFC_FFI_LC_SLOPES_DOWN          1
#define RFC_FFI_LC_SLOPES_ALL           2
#define RFC_FFI_LC_FVA                  3

/* Miner variants for damage_from_rp (rfc_rp_damage_method). */
#define RFC_FFI_RP_DEFAULT              0
#define RFC_FFI_RP_ELEMENTAR            1
#define RFC_FFI_RP_MODIFIED             2
#define RFC_FFI_RP_CONSISTENT           3

/* Counting methods (rfc_counting_method). */
#define RFC_FFI_COUNT_4PTM              1
#define RFC_FFI_COUNT_HCM               2
#define RFC_FFI_COUNT_ASTM              3

/* Flag bits. Same positions as rfc_flags. RFC_FFI_FLAGS_DEFAULT resolves inside init. */
#define RFC_FFI_FLAG_COUNT_RFM              (1 << 0)
#define RFC_FFI_FLAG_COUNT_DAMAGE           (1 << 1)
#define RFC_FFI_FLAG_COUNT_DH               (1 << 2)
#define RFC_FFI_FLAG_COUNT_RP               (1 << 3)
#define RFC_FFI_FLAG_COUNT_LC_UP            (1 << 4)
#define RFC_FFI_FLAG_COUNT_LC_DN            (1 << 5)
#define RFC_FFI_FLAG_COUNT_LC               ((1 << 4) | (1 << 5))
#define RFC_FFI_FLAG_COUNT_MK               (1 << 6)
#define RFC_FFI_FLAG_ENFORCE_MARGIN         (1 << 7)
#define RFC_FFI_FLAG_TPPRUNE_PRESERVE_POS   (1 << 8)
#define RFC_FFI_FLAG_TPPRUNE_PRESERVE_RES   (1 << 9)
#define RFC_FFI_FLAG_TPAUTOPRUNE            (1 << 10)
#define RFC_FFI_FLAG_AUTORESIZE             (1 << 11)
#define RFC_FFI_FLAGS_DEFAULT               (-1)

/* Counting state (rfc_state). */
#define RFC_FFI_STATE_INIT0             0
#define RFC_FFI_STATE_INIT              1
#define RFC_FFI_STATE_BUSY              2
#define RFC_FFI_STATE_BUSY_INTERIM      3
#define RFC_FFI_STATE_FINALIZE          4
#define RFC_FFI_STATE_FINISHED          5
#define RFC_FFI_STATE_ERROR             6

/* Wöhler parameters. Field order matches rfc_wl_param, natural alignment. */
typedef struct rfc_ffi_wl
{
    double sd;
    double nd;
    double k;
    double sx;
    double nx;
    double k2;
    double omission;
    double q;
    double q2;
    double D;
} rfc_ffi_wl_s;

/* One turning point. pos and adj_pos are 1-based sample positions. */
typedef struct rfc_ffi_tp
{
    double   value;
    double   damage;
    uint64_t pos;
    uint64_t adj_pos;
} rfc_ffi_tp_s;

/* Residue point before the pending-cycle strip. cls is the 0-based class. */
typedef struct rfc_ffi_residue
{
    double   value;
    uint32_t cls;
} rfc_ffi_residue_s;

/* Create a zeroed context and set its version word. NULL on allocation failure.
 * Destroy every handle, including one whose init failed. */
rfc_ffi *rfc_ffi_create(void);
void     rfc_ffi_destroy(rfc_ffi *h);

/* class_width must be > 0 when class_count > 0. flags may be RFC_FFI_FLAGS_DEFAULT.
 * Turning-point storage is allocated when the core was built with TP support. */
int32_t rfc_ffi_init(rfc_ffi *h, uint32_t class_count, double class_width,
                     double class_offset, double hysteresis, int32_t flags);
int32_t rfc_ffi_deinit(rfc_ffi *h);

int32_t rfc_ffi_set_flags(rfc_ffi *h, int32_t flags);
int32_t rfc_ffi_get_flags(const rfc_ffi *h, int32_t *flags);

/* method is RFC_FFI_LC_*. Adjusts the LC flag bits and the conversion mode.
 * Call while the context is still in RFC_FFI_STATE_INIT. */
int32_t rfc_ffi_set_lc_method(rfc_ffi *h, int32_t method);

/* method is RFC_FFI_COUNT_*. HCM and ASTM return UNSUPPORTED when compiled out. */
int32_t rfc_ffi_set_counting_method(rfc_ffi *h, int32_t method);

int32_t rfc_ffi_wl_init_modified(rfc_ffi *h, double sx, double nx, double k, double k2);
int32_t rfc_ffi_wl_init_any(rfc_ffi *h, const rfc_ffi_wl_s *wl);
int32_t rfc_ffi_wl_param_get(const rfc_ffi *h, rfc_ffi_wl_s *out);
int32_t rfc_ffi_wl_param_get_impaired(const rfc_ffi *h, rfc_ffi_wl_s *out);

/* Allocate a damage-history buffer of `cap` samples. method is RFC_FFI_SD_*. */
int32_t rfc_ffi_dh_init(rfc_ffi *h, int32_t method, uint64_t cap);

int32_t rfc_ffi_feed(rfc_ffi *h, const double *data, uint64_t count);
int32_t rfc_ffi_finalize(rfc_ffi *h, int32_t residual_method);

int32_t rfc_ffi_class_count(const rfc_ffi *h, uint32_t *count);
int32_t rfc_ffi_damage(const rfc_ffi *h, double *damage, double *damage_residue);

/* counts and amplitude each have class_count elements. amplitude is Sa, not range. */
int32_t rfc_ffi_rp(const rfc_ffi *h, double *counts, double *amplitude);
/* counts and level each have class_count elements. level is the class upper bound. */
int32_t rfc_ffi_lc(const rfc_ffi *h, double *counts, double *level);
/* Row-major class_count^2 raw counts (not divided by the full-cycle increment). */
int32_t rfc_ffi_rfm(const rfc_ffi *h, double *rfm_row_major);

/* Buffer protocol for variable-length results:
 *   buf == NULL: write the required element count to *count and return OK.
 *   *count < needed: write the required count and return INVARG.
 *   otherwise: fill buf and set *count to the number of elements written.
 */
int32_t rfc_ffi_residue(const rfc_ffi *h, rfc_ffi_residue_s *buf, uint32_t *count);
int32_t rfc_ffi_tp(const rfc_ffi *h, rfc_ffi_tp_s *buf, uint64_t *count);
int32_t rfc_ffi_dh(const rfc_ffi *h, double *buf, uint64_t *count);

int32_t     rfc_ffi_state(const rfc_ffi *h);
int32_t     rfc_ffi_error(const rfc_ffi *h);
const char *rfc_ffi_error_string(int32_t err);

/* Full-cycle count increment of this build (1.0, or 2.0 for integral counts). */
double rfc_ffi_full_cycle_increment(void);

/* Deep copy of the C counting state, including residue, histograms, turning
 * points, damage history, and the internal fields a preview finalize needs.
 * COUNT_DH is cleared on the copy so finalize does not require the original
 * sample stream. External amplitude-transform reference pointers are aliased
 * and must outlive the clone. On success *out is a new handle; destroy it.
 *
 * Python and Rust call this for a preview finalize. Turning points are the
 * C buffer inside the handle. */
int32_t rfc_ffi_clone(const rfc_ffi *src, rfc_ffi **out);

/* Sa and Sm are a Haigh reference curve, or both NULL with count 0 for FKM.
 * The pointers are stored and must stay valid for the life of the handle
 * (and of any clone) when they do not point at the internal curve. */
int32_t rfc_ffi_at_init(rfc_ffi *h, const double *Sa, const double *Sm, uint32_t count,
                        double M, double Sm_rig, double R_rig,
                        int32_t R_pinned, int32_t symmetric);
int32_t rfc_ffi_at_transform(const rfc_ffi *h, double Sa, double Sm, double *Sa_out);

/* counts and Sa have class_count elements and Sa is ascending. */
int32_t rfc_ffi_damage_from_rp(rfc_ffi *h, double *damage,
                               const double *counts, const double *Sa, int32_t method);

#ifdef __cplusplus
}
#endif

#endif /* RAINFLOW_FFI_H */
