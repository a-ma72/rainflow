/*
 * FFI implementation. Feature #ifdefs stay in this file.
 * The clone copies the C counting state, including the turning-point buffer.
 * Python and Rust both call rfc_ffi_clone() for a preview finalize.
 */

#include "rainflow_ffi.h"
#include "rainflow.h"

#include <stdlib.h>
#include <string.h>

typedef char rfc_ffi_require_double_value[(sizeof(rfc_value_t) == sizeof(double)) ? 1 : -1];

struct rfc_ffi
{
    rfc_ctx_s ctx;
};

static int32_t ffi_fail(rfc_ctx_s *ctx)
{
    rfc_error_e err;

    if (!ctx || ctx->version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }

    err = ctx->error;
    if (err == RFC_ERROR_NOERROR)
    {
        return RFC_FFI_ERR_UNEXP;
    }
    return (int32_t)err;
}

static int32_t ffi_ok(rfc_ctx_s *ctx, bool ok)
{
    return ok ? RFC_FFI_OK : ffi_fail(ctx);
}

static int32_t ffi_require_init_state(rfc_ffi *h)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (h->ctx.state != RFC_STATE_INIT)
    {
        return RFC_FFI_ERR_INVARG;
    }
    return RFC_FFI_OK;
}

#if !RFC_MINIMAL
static void ffi_wl_from_param(rfc_ffi_wl_s *out, const rfc_wl_param_s *wl)
{
    out->sd = wl->sd;
    out->nd = wl->nd;
    out->k = wl->k;
    out->sx = wl->sx;
    out->nx = wl->nx;
    out->k2 = wl->k2;
    out->omission = wl->omission;
    out->q = wl->q;
    out->q2 = wl->q2;
    out->D = wl->D;
}

static void ffi_wl_to_param(rfc_wl_param_s *wl, const rfc_ffi_wl_s *in)
{
    memset(wl, 0, sizeof(*wl));
    wl->sd = in->sd;
    wl->nd = in->nd;
    wl->k = in->k;
    wl->sx = in->sx;
    wl->nx = in->nx;
    wl->k2 = in->k2;
    wl->omission = in->omission;
    wl->q = in->q;
    wl->q2 = in->q2;
    wl->D = in->D;
}
#endif /* !RFC_MINIMAL */

rfc_ffi *rfc_ffi_create(void)
{
    rfc_ffi *h = (rfc_ffi *)calloc(1, sizeof(*h));
    if (!h)
    {
        return NULL;
    }
    h->ctx.version = sizeof(rfc_ctx_s);
    return h;
}

void rfc_ffi_destroy(rfc_ffi *h)
{
    if (!h)
    {
        return;
    }
    if (h->ctx.version == sizeof(rfc_ctx_s) && h->ctx.state >= RFC_STATE_INIT)
    {
        RFC_deinit(&h->ctx);
    }
    free(h);
}

int32_t rfc_ffi_deinit(rfc_ffi *h)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (h->ctx.state < RFC_STATE_INIT)
    {
        return RFC_FFI_OK;
    }
    return ffi_ok(&h->ctx, RFC_deinit(&h->ctx));
}

int32_t rfc_ffi_init(rfc_ffi *h, uint32_t class_count, double class_width,
                     double class_offset, double hysteresis, int32_t flags)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (class_count > RFC_CLASS_COUNT_MAX)
    {
        return RFC_FFI_ERR_INVARG;
    }

    if (!RFC_init(&h->ctx, (unsigned)class_count, class_width, class_offset,
                  hysteresis, (rfc_flags_e)flags))
    {
        return ffi_fail(&h->ctx);
    }

#if RFC_TP_SUPPORT
    if (!RFC_tp_init(&h->ctx, NULL, 64, false))
    {
        return ffi_fail(&h->ctx);
    }
#endif
    return RFC_FFI_OK;
}

int32_t rfc_ffi_set_flags(rfc_ffi *h, int32_t flags)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
#if RFC_MINIMAL
    (void)flags;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    return ffi_ok(&h->ctx, RFC_flags_set(&h->ctx, flags, RFC_FLAG_STACK_MAIN, true));
#endif
}

int32_t rfc_ffi_get_flags(const rfc_ffi *h, int32_t *flags)
{
    int got = 0;

    if (!h || !flags || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    (void)got;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    if (!RFC_flags_get(&h->ctx, &got, RFC_FLAG_STACK_MAIN))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    *flags = got;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_set_lc_method(rfc_ffi *h, int32_t method)
{
    int flags = 0;
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
    if (method < RFC_FFI_LC_SLOPES_UP || method > RFC_FFI_LC_FVA)
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    (void)flags;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    if (!RFC_flags_get(&h->ctx, &flags, RFC_FLAG_STACK_MAIN))
    {
        return ffi_fail(&h->ctx);
    }
    flags &= ~RFC_FLAGS_COUNT_LC;
    switch (method)
    {
        case RFC_FFI_LC_SLOPES_UP:
            flags |= RFC_FLAGS_COUNT_LC_UP;
            break;
        case RFC_FFI_LC_SLOPES_DOWN:
            flags |= RFC_FLAGS_COUNT_LC_DN;
            break;
        case RFC_FFI_LC_SLOPES_ALL:
        case RFC_FFI_LC_FVA:
            flags |= RFC_FLAGS_COUNT_LC;
            break;
        default:
            return RFC_FFI_ERR_INVARG;
    }
    if (!RFC_flags_set(&h->ctx, flags, RFC_FLAG_STACK_MAIN, true))
    {
        return ffi_fail(&h->ctx);
    }
    h->ctx.lc_count_method = (rfc_lc_count_method_e)method;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_set_counting_method(rfc_ffi *h, int32_t method)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
#if RFC_MINIMAL
    (void)method;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    if (method == RFC_FFI_COUNT_HCM)
    {
#if !RFC_HCM_SUPPORT
        return RFC_FFI_ERR_UNSUPPORTED;
#endif
    }
    else if (method == RFC_FFI_COUNT_ASTM)
    {
#if !RFC_ASTM_SUPPORT
        return RFC_FFI_ERR_UNSUPPORTED;
#endif
    }
    else if (method != RFC_FFI_COUNT_4PTM)
    {
        return RFC_FFI_ERR_INVARG;
    }
    h->ctx.counting_method = (rfc_counting_method_e)method;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_wl_init_modified(rfc_ffi *h, double sx, double nx, double k, double k2)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
#if RFC_MINIMAL
    (void)sx;
    (void)nx;
    (void)k;
    (void)k2;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    return ffi_ok(&h->ctx, RFC_wl_init_modified(&h->ctx, sx, nx, k, k2));
#endif
}

int32_t rfc_ffi_wl_init_any(rfc_ffi *h, const rfc_ffi_wl_s *wl)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
    if (!wl)
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    {
        rfc_wl_param_s param;
        ffi_wl_to_param(&param, wl);
        return ffi_ok(&h->ctx, RFC_wl_init_any(&h->ctx, &param));
    }
#endif
}

int32_t rfc_ffi_wl_param_get(const rfc_ffi *h, rfc_ffi_wl_s *out)
{
    if (!h || !out || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    {
        rfc_wl_param_s param;
        if (!RFC_wl_param_get(&h->ctx, &param))
        {
            return ffi_fail((rfc_ctx_s *)&h->ctx);
        }
        ffi_wl_from_param(out, &param);
        return RFC_FFI_OK;
    }
#endif
}

int32_t rfc_ffi_wl_param_get_impaired(const rfc_ffi *h, rfc_ffi_wl_s *out)
{
    if (!h || !out || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    {
        rfc_wl_param_s param;
        if (!RFC_wl_param_get_impaired(&h->ctx, &param))
        {
            return ffi_fail((rfc_ctx_s *)&h->ctx);
        }
        ffi_wl_from_param(out, &param);
        return RFC_FFI_OK;
    }
#endif
}

int32_t rfc_ffi_dh_init(rfc_ffi *h, int32_t method, uint64_t cap)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
    if (method < RFC_FFI_SD_NONE || method > RFC_FFI_SD_TRANSIENT_23c)
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (cap > (uint64_t)SIZE_MAX)
    {
        return RFC_FFI_ERR_INVARG;
    }
#if !RFC_DH_SUPPORT
    (void)method;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    return ffi_ok(&h->ctx, RFC_dh_init(&h->ctx, (rfc_sd_method_e)method, NULL, (size_t)cap, false));
#endif
}

int32_t rfc_ffi_feed(rfc_ffi *h, const double *data, uint64_t count)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (count > (uint64_t)SIZE_MAX)
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (!data && count != 0)
    {
        return RFC_FFI_ERR_INVARG;
    }
    return ffi_ok(&h->ctx, RFC_feed(&h->ctx, data, (size_t)count));
}

int32_t rfc_ffi_finalize(rfc_ffi *h, int32_t residual_method)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (residual_method < RFC_FFI_RES_NONE || residual_method > RFC_FFI_RES_RP_DIN45667)
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    if (residual_method > RFC_FFI_RES_NO_FINALIZE)
    {
        return RFC_FFI_ERR_UNSUPPORTED;
    }
#endif
    return ffi_ok(&h->ctx, RFC_finalize(&h->ctx, (rfc_res_method_e)residual_method));
}

int32_t rfc_ffi_class_count(const rfc_ffi *h, uint32_t *count)
{
    unsigned n = 0;

    if (!h || !count || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
#if RFC_MINIMAL
    if (h->ctx.state < RFC_STATE_INIT)
    {
        return RFC_FFI_ERR_INVARG;
    }
    *count = (uint32_t)h->ctx.class_count;
    return RFC_FFI_OK;
#else
    if (!RFC_class_count(&h->ctx, &n))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    *count = (uint32_t)n;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_damage(const rfc_ffi *h, double *damage, double *damage_residue)
{
    if (!h || !damage || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (h->ctx.state < RFC_STATE_INIT || h->ctx.state > RFC_STATE_FINISHED)
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
#if RFC_MINIMAL
    *damage = h->ctx.damage;
    if (damage_residue)
    {
        *damage_residue = h->ctx.damage_residue;
    }
    return RFC_FFI_OK;
#else
    return ffi_ok((rfc_ctx_s *)&h->ctx, RFC_damage(&h->ctx, damage, damage_residue));
#endif
}

#if !RFC_MINIMAL
static void ffi_copy_counts(double *dst, const rfc_counts_t *src, unsigned n)
{
    unsigned i;
    for (i = 0; i < n; i++)
    {
        dst[i] = (double)src[i];
    }
}
#endif

int32_t rfc_ffi_rp(const rfc_ffi *h, double *counts, double *amplitude)
{
#if RFC_MINIMAL
    (void)h;
    (void)counts;
    (void)amplitude;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    unsigned n = 0;
    rfc_counts_t *ct;
    rfc_value_t *sa;
    int32_t rc;

    if (!h || !counts || !amplitude || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (!RFC_class_count(&h->ctx, &n))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    if (n == 0)
    {
        return RFC_FFI_OK;
    }
    ct = (rfc_counts_t *)calloc(n, sizeof(*ct));
    sa = (rfc_value_t *)calloc(n, sizeof(*sa));
    if (!ct || !sa)
    {
        free(ct);
        free(sa);
        return RFC_FFI_ERR_MEMORY;
    }
    if (!RFC_rp_get(&h->ctx, ct, sa))
    {
        rc = ffi_fail((rfc_ctx_s *)&h->ctx);
        free(ct);
        free(sa);
        return rc;
    }
    ffi_copy_counts(counts, ct, n);
    ffi_copy_counts(amplitude, sa, n);
    free(ct);
    free(sa);
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_lc(const rfc_ffi *h, double *counts, double *level)
{
#if RFC_MINIMAL
    (void)h;
    (void)counts;
    (void)level;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    unsigned n = 0;
    rfc_counts_t *ct;
    rfc_value_t *lv;
    int32_t rc;

    if (!h || !counts || !level || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (!RFC_class_count(&h->ctx, &n))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    if (n == 0)
    {
        return RFC_FFI_OK;
    }
    ct = (rfc_counts_t *)calloc(n, sizeof(*ct));
    lv = (rfc_value_t *)calloc(n, sizeof(*lv));
    if (!ct || !lv)
    {
        free(ct);
        free(lv);
        return RFC_FFI_ERR_MEMORY;
    }
    if (!RFC_lc_get(&h->ctx, ct, lv))
    {
        rc = ffi_fail((rfc_ctx_s *)&h->ctx);
        free(ct);
        free(lv);
        return rc;
    }
    ffi_copy_counts(counts, ct, n);
    ffi_copy_counts(level, lv, n);
    free(ct);
    free(lv);
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_rfm(const rfc_ffi *h, double *rfm_row_major)
{
    unsigned n;
    size_t cells;
    size_t i;

    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (h->ctx.state < RFC_STATE_INIT)
    {
        return RFC_FFI_ERR_INVARG;
    }
    n = h->ctx.class_count;
    if (n == 0)
    {
        return RFC_FFI_OK;
    }
    if (!h->ctx.rfm)
    {
        return RFC_FFI_ERR_UNSUPPORTED;
    }
    if (!rfm_row_major)
    {
        return RFC_FFI_ERR_INVARG;
    }
    cells = (size_t)n * (size_t)n;
    for (i = 0; i < cells; i++)
    {
        rfm_row_major[i] = (double)h->ctx.rfm[i];
    }
    return RFC_FFI_OK;
}

int32_t rfc_ffi_residue(const rfc_ffi *h, rfc_ffi_residue_s *buf, uint32_t *count)
{
    const rfc_value_tuple_s *residue = NULL;
    unsigned n = 0;
    unsigned i;

    if (!h || !count || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (!RFC_res_get(&h->ctx, &residue, &n))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    if (!buf)
    {
        *count = n;
        return RFC_FFI_OK;
    }
    if (*count < n)
    {
        *count = n;
        return RFC_FFI_ERR_INVARG;
    }
    for (i = 0; i < n; i++)
    {
        buf[i].value = residue[i].value;
        buf[i].cls = (uint32_t)residue[i].cls;
    }
    *count = n;
    return RFC_FFI_OK;
}

int32_t rfc_ffi_tp(const rfc_ffi *h, rfc_ffi_tp_s *buf, uint64_t *count)
{
#if !RFC_TP_SUPPORT
    (void)h;
    (void)buf;
    if (count)
    {
        *count = 0;
    }
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    uint64_t n;
    uint64_t i;

    if (!h || !count || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (h->ctx.state < RFC_STATE_INIT)
    {
        return RFC_FFI_ERR_INVARG;
    }
    n = (uint64_t)h->ctx.tp_cnt;
    if (!buf)
    {
        *count = n;
        return RFC_FFI_OK;
    }
    if (*count < n)
    {
        *count = n;
        return RFC_FFI_ERR_INVARG;
    }
    for (i = 0; i < n; i++)
    {
        const rfc_value_tuple_s *tp = &h->ctx.tp[i];
        buf[i].value = tp->value;
        buf[i].pos = (uint64_t)tp->pos;
        buf[i].adj_pos = (uint64_t)tp->adj_pos;
#if RFC_DH_SUPPORT
        buf[i].damage = tp->damage;
#else
        buf[i].damage = 0.0;
#endif
    }
    *count = n;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_dh(const rfc_ffi *h, double *buf, uint64_t *count)
{
#if !RFC_DH_SUPPORT
    (void)h;
    (void)buf;
    if (count)
    {
        *count = 0;
    }
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    const double *dh = NULL;
    size_t n = 0;
    size_t i;

    if (!h || !count || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (!RFC_dh_get(&h->ctx, &dh, &n))
    {
        return ffi_fail((rfc_ctx_s *)&h->ctx);
    }
    if (!dh)
    {
        n = 0;
    }
    if (!buf)
    {
        *count = (uint64_t)n;
        return RFC_FFI_OK;
    }
    if (*count < (uint64_t)n)
    {
        *count = (uint64_t)n;
        return RFC_FFI_ERR_INVARG;
    }
    for (i = 0; i < n; i++)
    {
        buf[i] = dh[i];
    }
    *count = (uint64_t)n;
    return RFC_FFI_OK;
#endif
}

int32_t rfc_ffi_state(const rfc_ffi *h)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_UNEXP;
    }
    return (int32_t)RFC_state_get(&h->ctx);
}

int32_t rfc_ffi_error(const rfc_ffi *h)
{
    if (!h || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    return (int32_t)RFC_error_get(&h->ctx);
}

const char *rfc_ffi_error_string(int32_t err)
{
    switch (err)
    {
        case RFC_FFI_OK:                    return "No error";
        case RFC_FFI_ERR_INVARG:            return "Invalid arguments passed";
        case RFC_FFI_ERR_UNSUPPORTED:       return "Unsupported feature";
        case RFC_FFI_ERR_MEMORY:            return "Error on memory allocation";
        case RFC_FFI_ERR_TP:                return "Error while processing turning points";
        case RFC_FFI_ERR_AT:                return "Error while amplitude transformation";
        case RFC_FFI_ERR_DH_BAD_STREAM:     return "Input stream must be unique";
        case RFC_FFI_ERR_DH:                return "Error while damage history calculation/access";
        case RFC_FFI_ERR_LUT:               return "Error while accessing look up tables";
        case RFC_FFI_ERR_DATA_OUT_OF_RANGE: return "Input data leaves class range";
        case RFC_FFI_ERR_DATA_INCONSISTENT: return "Processed data is inconsistent";
        default:                            return "Unexpected error";
    }
}

double rfc_ffi_full_cycle_increment(void)
{
    return (double)RFC_FULL_CYCLE_INCREMENT;
}

int32_t rfc_ffi_at_init(rfc_ffi *h, const double *Sa, const double *Sm, uint32_t count,
                        double M, double Sm_rig, double R_rig,
                        int32_t R_pinned, int32_t symmetric)
{
    int32_t gate = ffi_require_init_state(h);
    if (gate != RFC_FFI_OK)
    {
        return gate;
    }
#if !RFC_AT_SUPPORT
    (void)Sa;
    (void)Sm;
    (void)count;
    (void)M;
    (void)Sm_rig;
    (void)R_rig;
    (void)R_pinned;
    (void)symmetric;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    return ffi_ok(&h->ctx, RFC_at_init(&h->ctx, Sa, Sm, (unsigned)count, M, Sm_rig, R_rig,
                                       R_pinned ? true : false, symmetric ? true : false));
#endif
}

int32_t rfc_ffi_at_transform(const rfc_ffi *h, double Sa, double Sm, double *Sa_out)
{
    if (!h || !Sa_out || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
#if !RFC_AT_SUPPORT
    (void)Sa;
    (void)Sm;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    return ffi_ok((rfc_ctx_s *)&h->ctx, RFC_at_transform(&h->ctx, Sa, Sm, Sa_out));
#endif
}

int32_t rfc_ffi_damage_from_rp(rfc_ffi *h, double *damage,
                               const double *counts, const double *Sa, int32_t method)
{
#if RFC_MINIMAL
    (void)h;
    (void)damage;
    (void)counts;
    (void)Sa;
    (void)method;
    return RFC_FFI_ERR_UNSUPPORTED;
#else
    unsigned n;
    unsigned i;
    rfc_counts_t *ct = NULL;
    rfc_value_t *sa = NULL;
    int32_t rc;

    if (!h || !damage || h->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (method < RFC_FFI_RP_DEFAULT || method > RFC_FFI_RP_CONSISTENT)
    {
        return RFC_FFI_ERR_INVARG;
    }
    n = h->ctx.class_count;
    if (!counts || !Sa || n == 0)
    {
        return RFC_FFI_ERR_INVARG;
    }
    ct = (rfc_counts_t *)calloc(n, sizeof(*ct));
    sa = (rfc_value_t *)calloc(n, sizeof(*sa));
    if (!ct || !sa)
    {
        free(ct);
        free(sa);
        return RFC_FFI_ERR_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        if (counts[i] < 0.0)
        {
            free(ct);
            free(sa);
            return RFC_FFI_ERR_INVARG;
        }
        ct[i] = (rfc_counts_t)counts[i];
        sa[i] = Sa[i];
    }
    if (!RFC_damage_from_rp(&h->ctx, damage, ct, sa, (rfc_rp_damage_method_e)method))
    {
        rc = ffi_fail(&h->ctx);
        free(ct);
        free(sa);
        return rc;
    }
    free(ct);
    free(sa);
    return RFC_FFI_OK;
#endif
}

#if RFC_TP_SUPPORT
static int32_t ffi_copy_tp(rfc_ffi *dst, const rfc_ffi *src)
{
    size_t cap = src->ctx.tp_cap;
    rfc_value_tuple_s *grown;

    if (cap < src->ctx.tp_cnt)
    {
        cap = src->ctx.tp_cnt;
    }
    if (cap < 64)
    {
        cap = 64;
    }
    if (!RFC_tp_init(&dst->ctx, NULL, cap, false))
    {
        return ffi_fail(&dst->ctx);
    }
    if (src->ctx.tp && src->ctx.tp_cnt)
    {
        if (src->ctx.tp_cnt > dst->ctx.tp_cap)
        {
            grown = (rfc_value_tuple_s *)dst->ctx.mem_alloc(
                dst->ctx.tp, src->ctx.tp_cnt, sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TP);
            if (!grown)
            {
                return RFC_FFI_ERR_MEMORY;
            }
            dst->ctx.tp = grown;
            dst->ctx.tp_cap = src->ctx.tp_cnt;
        }
        memcpy(dst->ctx.tp, src->ctx.tp, src->ctx.tp_cnt * sizeof(rfc_value_tuple_s));
    }
    dst->ctx.tp_cnt = src->ctx.tp_cnt;
    dst->ctx.tp_locked = src->ctx.tp_locked;
    dst->ctx.tp_prune_size = src->ctx.tp_prune_size;
    dst->ctx.tp_prune_threshold = src->ctx.tp_prune_threshold;
    dst->ctx.internal.margin[0] = src->ctx.internal.margin[0];
    dst->ctx.internal.margin[1] = src->ctx.internal.margin[1];
    dst->ctx.internal.margin_stage = src->ctx.internal.margin_stage;
    return RFC_FFI_OK;
}
#endif

int32_t rfc_ffi_clone(const rfc_ffi *src, rfc_ffi **out)
{
    rfc_ffi *dst;
    int flags;
    size_t res_n;
    int32_t rc;
#if !RFC_MINIMAL
    rfc_wl_param_s wl;
#endif

    if (!src || !out || src->ctx.version != sizeof(rfc_ctx_s))
    {
        return RFC_FFI_ERR_INVARG;
    }
    if (src->ctx.state < RFC_STATE_INIT || src->ctx.state > RFC_STATE_FINISHED)
    {
        return RFC_FFI_ERR_INVARG;
    }

    dst = rfc_ffi_create();
    if (!dst)
    {
        return RFC_FFI_ERR_MEMORY;
    }

    flags = src->ctx.internal.flags;
#if RFC_DH_SUPPORT
    flags &= ~RFC_FLAGS_COUNT_DH;
#endif

    if (!RFC_init(&dst->ctx, src->ctx.class_count, src->ctx.class_width,
                  src->ctx.class_offset, src->ctx.hysteresis, (rfc_flags_e)flags))
    {
        rc = ffi_fail(&dst->ctx);
        rfc_ffi_destroy(dst);
        return rc;
    }

#if !RFC_MINIMAL
    if (!RFC_wl_param_get(&src->ctx, &wl) || !RFC_wl_init_any(&dst->ctx, &wl))
    {
        rc = ffi_fail(&dst->ctx);
        rfc_ffi_destroy(dst);
        return rc;
    }
    dst->ctx.counting_method = src->ctx.counting_method;
    dst->ctx.lc_count_method = src->ctx.lc_count_method;
    dst->ctx.series_start = src->ctx.series_start;
    dst->ctx.series_end = src->ctx.series_end;
    dst->ctx.series_bounds_valid = src->ctx.series_bounds_valid;
    dst->ctx.internal.wl = src->ctx.internal.wl;
#endif

    dst->ctx.full_inc = src->ctx.full_inc;
    dst->ctx.half_inc = src->ctx.half_inc;
    dst->ctx.curr_inc = src->ctx.curr_inc;
    dst->ctx.damage = src->ctx.damage;
    dst->ctx.damage_residue = src->ctx.damage_residue;

    if (src->ctx.rfm && dst->ctx.rfm)
    {
        memcpy(dst->ctx.rfm, src->ctx.rfm,
               (size_t)src->ctx.class_count * (size_t)src->ctx.class_count * sizeof(*dst->ctx.rfm));
    }
#if !RFC_MINIMAL
    if (src->ctx.rp && dst->ctx.rp)
    {
        memcpy(dst->ctx.rp, src->ctx.rp, (size_t)src->ctx.class_count * sizeof(*dst->ctx.rp));
    }
    if (src->ctx.lc && dst->ctx.lc)
    {
        memcpy(dst->ctx.lc, src->ctx.lc, (size_t)src->ctx.class_count * sizeof(*dst->ctx.lc));
    }
#endif

    res_n = src->ctx.residue_cnt + ((src->ctx.state == RFC_STATE_BUSY_INTERIM) ? 1u : 0u);
    if (res_n > dst->ctx.residue_cap)
    {
        rfc_ffi_destroy(dst);
        return RFC_FFI_ERR_MEMORY;
    }
    if (res_n && src->ctx.residue && dst->ctx.residue)
    {
        memcpy(dst->ctx.residue, src->ctx.residue, res_n * sizeof(*dst->ctx.residue));
    }
    dst->ctx.residue_cnt = src->ctx.residue_cnt;

    dst->ctx.internal.slope = src->ctx.internal.slope;
    dst->ctx.internal.extrema[0] = src->ctx.internal.extrema[0];
    dst->ctx.internal.extrema[1] = src->ctx.internal.extrema[1];
#if RFC_GLOBAL_EXTREMA
    dst->ctx.internal.extrema_changed = src->ctx.internal.extrema_changed;
#endif
    dst->ctx.internal.pos = src->ctx.internal.pos;
    dst->ctx.internal.pos_offset = src->ctx.internal.pos_offset;

#if RFC_TP_SUPPORT
    rc = ffi_copy_tp(dst, src);
    if (rc != RFC_FFI_OK)
    {
        rfc_ffi_destroy(dst);
        return rc;
    }
#endif

#if RFC_DH_SUPPORT
    if (src->ctx.dh)
    {
        size_t cap = src->ctx.dh_cap ? src->ctx.dh_cap : src->ctx.dh_cnt;
        if (cap)
        {
            if (!RFC_dh_init(&dst->ctx, src->ctx.spread_damage_method, NULL, cap, false))
            {
                rc = ffi_fail(&dst->ctx);
                rfc_ffi_destroy(dst);
                return rc;
            }
            if (src->ctx.dh_cnt)
            {
                memcpy(dst->ctx.dh, src->ctx.dh, src->ctx.dh_cnt * sizeof(double));
            }
            dst->ctx.dh_cnt = src->ctx.dh_cnt;
        }
    }
#endif

#if RFC_AT_SUPPORT
    dst->ctx.internal.at_haigh = src->ctx.internal.at_haigh;
    dst->ctx.at = src->ctx.at;
    if (src->ctx.at.Sa == src->ctx.internal.at_haigh.Sa)
    {
        dst->ctx.at.Sa = dst->ctx.internal.at_haigh.Sa;
        dst->ctx.at.Sm = dst->ctx.internal.at_haigh.Sm;
    }
#if RFC_DAMAGE_FAST
    /* dst->init rebuilt the LUT without the amplitude transform. Force the
     * slow path until a later at_init rebuilds the tables. */
    dst->ctx.damage_lut_inapt = 1;
#endif
#endif

#if RFC_HCM_SUPPORT
    dst->ctx.internal.hcm.IR = src->ctx.internal.hcm.IR;
    dst->ctx.internal.hcm.IZ = src->ctx.internal.hcm.IZ;
    if (src->ctx.internal.hcm.stack && src->ctx.internal.hcm.stack_cap)
    {
        if (!dst->ctx.internal.hcm.stack ||
            dst->ctx.internal.hcm.stack_cap < src->ctx.internal.hcm.stack_cap)
        {
            rfc_value_tuple_s *grown = (rfc_value_tuple_s *)dst->ctx.mem_alloc(
                dst->ctx.internal.hcm.stack, src->ctx.internal.hcm.stack_cap,
                sizeof(rfc_value_tuple_s), RFC_MEM_AIM_HCM);
            if (!grown)
            {
                rfc_ffi_destroy(dst);
                return RFC_FFI_ERR_MEMORY;
            }
            dst->ctx.internal.hcm.stack = grown;
            dst->ctx.internal.hcm.stack_cap = src->ctx.internal.hcm.stack_cap;
        }
        memcpy(dst->ctx.internal.hcm.stack, src->ctx.internal.hcm.stack,
               src->ctx.internal.hcm.stack_cap * sizeof(*dst->ctx.internal.hcm.stack));
    }
#endif

    dst->ctx.state = src->ctx.state;
    dst->ctx.error = src->ctx.error;
    *out = dst;
    return RFC_FFI_OK;
}
