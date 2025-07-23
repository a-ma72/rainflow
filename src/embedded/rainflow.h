/**
 * @file rainflow.h
 * @brief Public API of the generic fixed-point rainflow engine
 *        (4-point method) for targets without an FPU.
 *        Application-independent — class count, class bounds and the
 *        damage LUT come from rainflow_config.h.
 *
 * Two-stage processing per incoming sample:
 *
 *   1. Q28.4 domain: hysteresis / turning-point filter on the raw
 *      sample. A turning point is confirmed only after it exceeds the
 *      configured hysteresis threshold. This suppresses noise BEFORE
 *      quantization.
 *   2. Integer domain: a confirmed turning point is classified
 *      (Q28.4 -> RF_Class_t, RF_NUM_CLASSES bins) and pushed onto the
 *      residue stack. The 4-point algorithm itself then uses only
 *      integers (class indices) — no fixed-point arithmetic on the
 *      hot counting path.
 *
 * Deliberately no full time-series buffer is kept — only the
 * classified residue stack, which is also the only working buffer
 * required. The module is therefore suitable for a never-ending
 * input stream.
 *
 * Further simplifications versus full ASTM E1049 / HCM:
 *  - For classified values the maximum residue-stack depth is bounded
 *    by the algorithm itself to 2 * RF_NUM_CLASSES — overflow is
 *    excluded when RF_NUM_CLASSES is correct.
 *  - Individual cycles (from/to) are not emitted — only the resulting
 *    damage increment matters. Amplitude -> damage is a table lookup
 *    (RF_DAMAGE_LUT) at cycle closure, with no intermediate storage.
 *
 * Reference: DIN 45667 / FVA 4-point method (ASTM E 1049-85 is the
 * 3-point algorithm). Conceptually validated against a-ma72/rainflow
 * (github.com/a-ma72/rainflow).
 */
#ifndef RAINFLOW_H
#define RAINFLOW_H

#include <stdint.h>
#include <stdbool.h>

/* ------------------------------------------------------------------ */
/* Rainflow configuration (provides RF_Value_t, RF_NUM_CLASSES,
 * RF_CLASS_MIN/MAX, RF_DAMAGE_LUT — must precede all type defs).     */
#include "rainflow_config.h"

/* ------------------------------------------------------------------ */
/* Classification (stage 2: confirmed turning points only)            */
/* ------------------------------------------------------------------ */

/**
 * Class index of a classified turning point, 0 .. RF_NUM_CLASSES-1.
 * uint8_t is enough for RF_NUM_CLASSES <= 256; switch to uint16_t
 * for a larger class count.
 */
typedef uint8_t RF_Class_t;

/* ------------------------------------------------------------------ */
/* Residue-stack capacity (= only working buffer)                     */
/* ------------------------------------------------------------------ */

/**
 * Algorithmic upper bound on residue depth from the 4-point condition
 * on classified values. No safety margin needed: each class can appear
 * in the residue at most twice as an "open flank".
 *
 * HOLDS ONLY under the precondition hysteresis >= RF_CLASS_WIDTH
 * (see there and RF_Init()) — otherwise the residue alternation
 * invariant can be broken and this bound is no longer reliable.
 */
#define RF_MAX_RESIDUUM (2U * RF_NUM_CLASSES)

/* ------------------------------------------------------------------ */
/* Status codes                                                       */
/* ------------------------------------------------------------------ */

typedef enum
{
    RF_OK                  = 0,
    RF_ERR_NULL_PTR        = -1,
    RF_ERR_RESIDUUM_FULL   = -2, /* see RF_ERR_INVALID_CONFIG — when
                                     hysteresis >= RF_CLASS_WIDTH is
                                     correctly enforced this should
                                     not occur; kept as a safety net
                                     against other config errors. */
    RF_ERR_INVALID_CONFIG  = -3  /* from RF_Init() when
                                     hysteresis < RF_CLASS_WIDTH —
                                     see RF_Init() for the rationale. */
} RF_Status_t;

/* ------------------------------------------------------------------ */
/* Direction of the Q28.4 peak/valley candidate                       */
/* ------------------------------------------------------------------ */

typedef enum
{
    RF_SLOPE_UNKNOWN = 0, /* no confirmed turning point yet */
    RF_SLOPE_RISING  = 1, /* running_extreme is a peak candidate */
    RF_SLOPE_FALLING = 2  /* running_extreme is a valley candidate */
} RF_Slope_t;

/* ------------------------------------------------------------------ */
/* 96-bit damage accumulator (12 bytes of value)                      */
/* ------------------------------------------------------------------ */

/**
 * Unsigned 96-bit value for accumulated total damage, composed of
 * 8 bytes (lo, uint64_t) and the remaining 4 bytes (hi, uint32_t) —
 * C has no native 96-bit integer. Range: 0 .. (2^96 - 1) = lo + hi * 2^64.
 *
 * @note sizeof(RF_Damage96_t) may be larger than 12 bytes due to
 *       uint64_t alignment (typically 16 bytes with 4 padding bytes
 *       after hi) — the 12 bytes refer to value-range bit depth, not
 *       the physical size of the struct.
 */
typedef struct
{
    uint64_t lo; /* bits 0..63 */
    uint32_t hi; /* bits 64..95 */
} RF_Damage96_t;

/* ------------------------------------------------------------------ */
/* Main context (fully static, no heap)                               */
/* ------------------------------------------------------------------ */

/**
 * @note NOT thread-safe. The context must never be accessed in
 *       parallel from multiple tasks/threads — all API functions are
 *       non-reentrant and require exclusive access.
 */
typedef struct
{
    /* --- Stage 1: hysteresis / turning-point filter, Q28.4 --- */
    RF_Value_t hysteresis;         /* threshold in Q28.4 */
    RF_Value_t running_extreme;    /* current peak/valley candidate */
    bool       has_running_extreme;
    RF_Slope_t slope;              /* direction of running_extreme */

    /* --- Stage 2: classified 4-point stack, integer domain --- */
    RF_Class_t residuum[RF_MAX_RESIDUUM];
    uint16_t   residuum_count;

    /* --- Accumulated total damage --- */
    /* Sum of all actually applied (not merely predicted) damage
     * increments: raised by RF_ProcessSample() on every closed cycle
     * and by RF_FlushResiduumRepeated() when commit == true. When
     * commit == false (Predict) it is left unchanged, because by
     * definition nothing is actually applied. 96-bit range (see
     * RF_Damage96_t) — practically never overflows for a single
     * increment (max uint32_t per call). Read via
     * RF_GetAccumulatedDamage().
     */
    RF_Damage96_t accumulated_damage;

    /* --- Range-pair counts (histogram over cycle range) --- */
    /* rp_counts[i] = number of actually closed cycles with range i
     * (i = class difference of the two closing points,
     * 0 .. RF_NUM_CLASSES-1 — the same index used for RF_DAMAGE_LUT).
     * Raised by RF_ProcessSample() on every closed cycle and by
     * RF_FlushResiduumRepeated() when commit == true. When
     * commit == false (Predict) it is left unchanged, same as
     * accumulated_damage. May theoretically overflow (uint32_t per
     * class) over a very long run — monitoring is the caller's
     * responsibility (see RF_GetRangePairCounts()).
     * Due to the hysteresis precondition (hysteresis >= RF_CLASS_WIDTH)
     * the first entry is guaranteed to stay 0.
     * The field is not required in production target code and may be
     * omitted to save RAM.
     */
    uint32_t rp_counts[RF_NUM_CLASSES];  /*optional*/
} RF_Ctx_t;

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

/**
 * @brief Initialize the rainflow context.
 *
 * @param[out] ctx        Pointer to the context to initialize.
 * @param[in]  hysteresis Hysteresis threshold in Q28.4 (see
 *                         RF_FIXED_SHIFT), applied to the raw sample
 *                         BEFORE classification.
 *
 *                         PRECONDITION: hysteresis >= RF_CLASS_WIDTH.
 *                         Only then are two consecutive confirmed
 *                         turning points (stage 1) guaranteed not to
 *                         fall into the same class after classification
 *                         (stage 2). If hysteresis is smaller, a
 *                         confirmed turning point may classify into the
 *                         same class as its predecessor; the plateau
 *                         guard in rf_stack_push_and_close_buf() then
 *                         discards that point entirely and with it the
 *                         information that a real reversal occurred —
 *                         the next confirmed point can then point in
 *                         the same direction relative to the (older)
 *                         stack top instead of the expected opposite
 *                         direction. That breaks not only the strict
 *                         residue alternation assumed in several
 *                         function comments, but in particular the
 *                         correctness proof of the wrap-around handling
 *                         in RF_FlushResiduumRepeated() (which
 *                         explicitly assumes a plateau can occur ONLY
 *                         at the wrap). RF_Init() therefore enforces
 *                         the precondition rather than merely documenting
 *                         it.
 *
 * @return RF_OK on success, RF_ERR_NULL_PTR if ctx is NULL,
 *         RF_ERR_INVALID_CONFIG if hysteresis < RF_CLASS_WIDTH
 *         (ctx is still set to a defined zeroed state — see rainflow.c —
 *         but must not be used because of the error code).
 */
RF_Status_t RF_Init(RF_Ctx_t *ctx, RF_Value_t hysteresis);

/**
 * @brief Process a new sample from the stream.
 *
 * Flow: (1) The Q28.4 hysteresis filter checks whether a turning point
 * is confirmed; if not, the function returns with
 * *out_damage_increment = 0. (2) A confirmed turning point is
 * classified and pushed onto the residue stack; the 4-point closure
 * condition is checked in cascade. (3) For each cycle closed, the
 * damage contribution from RF_DAMAGE_LUT is added AND the counter for
 * that class difference in ctx->rp_counts is incremented (see
 * RF_GetRangePairCounts()) — *out_damage_increment returns only the
 * summed damage of this call, not the individual class differences;
 * query RF_GetRangePairCounts() after the stream for the full
 * distribution. (4) This sum is ALWAYS also added to
 * ctx->accumulated_damage (see RF_GetAccumulatedDamage()) —
 * regardless of whether out_damage_increment is set.
 *
 * @param[in,out] ctx                  Rainflow context.
 * @param[in]     sample               New sample in Q28.4 fixed-point.
 * @param[out]    out_damage_increment Optional (may be NULL) if this
 *                                     call's damage increment is not
 *                                     needed separately — the total
 *                                     remains available via
 *                                     RF_GetAccumulatedDamage(). If
 *                                     given: sum of damage contributions
 *                                     of all cycles closed in this call
 *                                     (0 if none closed or no turning
 *                                     point was confirmed).
 *
 * @return RF_OK on success.
 *         RF_ERR_NULL_PTR if ctx is NULL.
 *         RF_ERR_RESIDUUM_FULL on a configuration error
 *         (see RF_MAX_RESIDUUM) — unreachable in normal operation.
 *
 * @pre RF_Init() must have been called successfully first.
 */
RF_Status_t RF_ProcessSample(RF_Ctx_t *ctx, RF_Value_t sample,
                              uint32_t *out_damage_increment);

/**
 * @brief Compute the damage increment that would result from closing
 *        the currently open residue remainder with the ASTM E1049
 *        "repeated residue" technique: the residue is notionally
 *        appended to itself (as if this load segment repeated
 *        periodically) and the 4-point algorithm is run once over
 *        that. This finds extra, larger/nested full cycles that would
 *        arise under periodic repetition of this load segment.
 *
 * Stage-1 candidate (running_extreme): the current, not-yet-confirmed
 * peak/valley candidate of hysteresis stage 1 (see RF_ProcessSample())
 * is included in the calculation — if it represents a point different
 * from the last confirmation — as if the signal had stopped right now.
 * Otherwise the last, possibly already long, branch would be ignored
 * systematically. Details and proof: see the comment in rainflow.c.
 *
 * @p commit controls whether the result is actually applied:
 *   - commit == true:  a "real" flush — afterwards the residue is
 *     reduced to the last point of the (possibly running_extreme-
 *     extended) sequence, kept as the base for continuing the stream.
 *     If running_extreme was included, it is treated as synthetically
 *     confirmed: stage-1 state is reset to "no candidate, no known
 *     direction" (same as the first-sample anchor case) —
 *     running_extreme itself is unchanged because its value already
 *     matches that new anchor.
 *   - commit == false: prediction only ("Predict") — the context is
 *     left completely unchanged, including residue, residuum_count,
 *     stage-1 state (running_extreme, slope) AND
 *     ctx->accumulated_damage. Useful e.g. before a real flush or
 *     purely informatively, without affecting the live counting state.
 *
 * When commit == true the computed damage increment is also added to
 * ctx->accumulated_damage (see RF_GetAccumulatedDamage()) and for each
 * cycle closed ctx->rp_counts[range] is incremented by one (see
 * RF_GetRangePairCounts()) — both regardless of whether
 * out_damage_increment is set; when commit == false both are omitted
 * as well.
 *
 * The Q28.4 raw value of running_extreme itself is never changed in
 * either case — only its role in stage-1 state (slope) can change when
 * commit == true (see above).
 *
 * Because the sample stream never really ends, this function must be
 * called periodically or on demand with commit == true (e.g. every N
 * seconds) so open remainder cycles are not left uncounted forever.
 *
 * @param[in,out] ctx                  Rainflow context. When
 *                                     commit == false used read-only
 *                                     (logically const).
 * @param[in]     commit               true: update residue (and
 *                                     possibly stage-1 state, see
 *                                     above) after the calculation.
 *                                     false: compute only, leave
 *                                     context completely unchanged
 *                                     (Predict).
 * @param[out]    out_damage_increment Optional (may be NULL) if this
 *                                     call's damage increment is not
 *                                     needed separately — when
 *                                     commit == true the total remains
 *                                     available via
 *                                     RF_GetAccumulatedDamage() (when
 *                                     commit == false, Predict, the
 *                                     value is then not available at
 *                                     all — that call only makes sense
 *                                     with out_damage_increment set).
 *                                     If given: sum of damage
 *                                     contributions of all (possibly
 *                                     only hypothetical) closed cycles.
 *
 * @return RF_OK on success, RF_ERR_NULL_PTR if ctx is NULL,
 *         RF_ERR_RESIDUUM_FULL on a configuration error — unreachable
 *         in normal operation.
 */
RF_Status_t RF_FlushResiduumRepeated(RF_Ctx_t *ctx, bool commit,
                                      uint32_t *out_damage_increment);

/**
 * @brief Return the current number of open points in the residue.
 *        Useful for diagnostics/monitoring.
 *
 * @param[in] ctx Rainflow context.
 * @return Number of open points, 0 if ctx is NULL.
 */
uint16_t RF_GetResiduumCount(const RF_Ctx_t *ctx);

/**
 * @brief Return accumulated total damage so far (96 bit).
 *
 * Sum of all actually applied (not merely predicted) damage increments
 * since the last RF_Init(): raised by RF_ProcessSample() on every
 * closed cycle and by RF_FlushResiduumRepeated() when commit == true.
 * RF_FlushResiduumRepeated() calls with commit == false (Predict) do
 * NOT contribute.
 *
 * @param[in] ctx Rainflow context.
 * @return Accumulated total damage as RF_Damage96_t
 *         ({lo=0, hi=0} if ctx is NULL).
 */
RF_Damage96_t RF_GetAccumulatedDamage(const RF_Ctx_t *ctx);

/**
 * @brief Return accumulated total damage so far as double.
 *
 * Convenience for further processing on targets that already use
 * floating point (e.g. logging/diagnostics on a PC tool, remaining-
 * life calculation). NOT intended for target code without an FPU —
 * see the module comment above: using this function on a target
 * without a hardware FPU pulls in expensive software float emulation
 * and contradicts the design of the module.
 *
 * Precision note: a double has only 53 bits of mantissa. Above 2^53
 * (~9.007e15) the low bits of lo are discarded when rounding to the
 * nearest representable double. For the exact 96-bit value see
 * RF_GetAccumulatedDamage().
 *
 * This is NOT a purely theoretical edge case: with RF_NUM_CLASSES=128
 * and S-N slope k=5 the largest possible damage increment of a single
 * classification (127^5) is already ~3.3e10. At 50 Hz sampling and
 * 8000 h of operation a conservative continuous-duty scenario
 * (order 1e16–1e19, depending on the share of fully closed maximum
 * cycles) can reach and exceed 2^53 — so the precision loss versus
 * RF_GetAccumulatedDamage() is realistic in long-running duty with
 * correspondingly sized S-N curves / damage LUTs, not only for
 * pathological extremes.
 *
 * @param[in] ctx Rainflow context.
 * @return Accumulated total damage as double (possibly rounded to
 *         double precision, see above), 0.0 if ctx is NULL.
 */
double RF_GetAccumulatedDamageDouble(const RF_Ctx_t *ctx);

/**
 * @brief Return the range-pair counts (rainflow histogram): number of
 *        actually closed cycles per class difference.
 *
 * Sum of all actually applied (not merely predicted) cycle closures
 * since the last RF_Init(), broken down by class difference — same
 * commit semantics as RF_GetAccumulatedDamage() (see there):
 * RF_ProcessSample() increments on every closed cycle,
 * RF_FlushResiduumRepeated() only when commit == true.
 *
 * @param[in] ctx Rainflow context.
 * @return Pointer to an array of RF_NUM_CLASSES uint32_t entries
 *         (index 0 .. RF_NUM_CLASSES-1 = class difference), valid
 *         as long as ctx lives and is not RF_Init()'d again.
 *         NULL if ctx is NULL.
 */
const uint32_t *RF_GetRangePairCounts(const RF_Ctx_t *ctx);

#endif /* RAINFLOW_H */
