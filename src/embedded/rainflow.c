/**
 * @file rainflow.c
 * @brief Implementation of the generic rainflow engine (see rainflow.h).
 *
 * The algorithm itself is application-independent and is not regenerated
 * per project; only rainflow_config.h changes with class grid and
 * S-N curve.
 *
 * Two stages:
 *  - Stage 1 (Q28.4): hysteresis / turning-point filter, see
 *    rf_hysteresis_filter().
 *  - Stage 2 (integer): classify confirmed turning points and close
 *    cycles on a 4-point stack with a direct damage-LUT lookup, see
 *    rf_classify() and rf_stack_push_and_close().
 */
#include "rainflow.h"
#include <string.h>

/* ------------------------------------------------------------------ */
/* Stage 1: Q28.4 hysteresis / turning-point filter                   */
/* ------------------------------------------------------------------ */

/**
 * @brief Process a raw sample in the hysteresis filter and report
 *        whether a turning point was confirmed.
 *
 * State machine with exactly one running peak/valley candidate
 * (ctx->running_extreme). The very first raw sample is always
 * confirmed immediately as the sequence start (required anchor for
 * stage 2). After that a new turning point is confirmed only when
 * the value moves by more than ctx->hysteresis against the current
 * candidate direction (strict, like rfcnt `delta > hysteresis`).
 *
 * Important invariant for stage 2 (see RF_FlushResiduumRepeated()):
 * because a new point is confirmed only after a move > ctx->hysteresis
 * against the candidate direction, two consecutive confirmed raw
 * values always differ strictly (no plateau in the raw domain). After
 * classification (rf_classify(), monotonically non-decreasing) rank
 * order is preserved, but a true reversal can never "flip" into a
 * continuation of the same direction through quantization — a plateau
 * (same class value) cannot occur internally in the residue.
 *
 * @param[in,out] ctx              Rainflow context (only stage-1 fields
 *                                 are read/written).
 * @param[in]     sample           New raw sample, Q28.4.
 * @param[out]    out_confirmed_pt On true return: the confirmed turning
 *                                 point (Q28.4). Unchanged on false.
 *
 * @return true if a turning point was confirmed in this call, else false.
 *
 * @pre ctx->hysteresis >= 0.
 */
static bool rf_hysteresis_filter(RF_Ctx_t *ctx, RF_Value_t sample,
                                  RF_Value_t *out_confirmed_pt)
{
    if (!ctx->has_running_extreme)
    {
        /* Very first raw sample: confirmed immediately as start anchor. */
        ctx->running_extreme = sample;
        ctx->has_running_extreme = true;
        ctx->slope = RF_SLOPE_UNKNOWN;
        *out_confirmed_pt = sample;
        return true;
    }

    const int32_t delta = (int32_t)sample - (int32_t)ctx->running_extreme;
    const int32_t hysteresis = (int32_t)ctx->hysteresis;

    if (ctx->slope == RF_SLOPE_UNKNOWN)
    {
        /* Direction not yet established: wait until the hysteresis
           threshold is exceeded in either direction. */
        /* rfcnt: delta > hysteresis (strict). Equality stays noise. */
        if ((delta <= hysteresis) && (delta >= -hysteresis))
        {
            return false;
        }

        ctx->slope = (delta > 0) ? RF_SLOPE_RISING : RF_SLOPE_FALLING;
        ctx->running_extreme = sample;
        return false;
    }

    if (ctx->slope == RF_SLOPE_RISING)
    {
        if (sample >= ctx->running_extreme)
        {
            ctx->running_extreme = sample; /* extend peak candidate */
            return false;
        }

        if (-delta <= hysteresis)
        {
            return false; /* reversal still inside the hysteresis band */
        }

        /* Reversal > hysteresis: peak is confirmed (like rfcnt). */
        *out_confirmed_pt = ctx->running_extreme;
        ctx->slope = RF_SLOPE_FALLING;
        ctx->running_extreme = sample;
        return true;
    }

    /* ctx->slope == RF_SLOPE_FALLING */
    if (sample <= ctx->running_extreme)
    {
        ctx->running_extreme = sample; /* extend valley candidate */
        return false;
    }

    if (delta <= hysteresis)
    {
        return false; /* reversal still inside the hysteresis band */
    }

    /* Reversal > hysteresis: valley is confirmed (like rfcnt). */
    *out_confirmed_pt = ctx->running_extreme;
    ctx->slope = RF_SLOPE_RISING;
    ctx->running_extreme = sample;
    return true;
}

/* ------------------------------------------------------------------ */
/* Stage 2: classification + 4-point stack, integer domain            */
/* ------------------------------------------------------------------ */

/**
 * @brief Classify a confirmed turning point (Q28.4) into a class
 *        index 0 .. RF_NUM_CLASSES-1.
 *
 * Floor binning like rfcnt `QUANTIZE`:
 *   class = (unsigned)((value - class_offset) / class_width)
 * in Q28.4: integer division without +width/2. Values below
 * RF_CLASS_MIN_FIXED fall into class 0, values at or above
 * RF_CLASS_MAX_FIXED into the last class (rfcnt would abort when
 * value >= offset + count*width).
 *
 * @param[in] value Raw value in Q28.4.
 * @return Class index 0 .. RF_NUM_CLASSES-1.
 */
static RF_Class_t rf_classify(RF_Value_t value)
{
    if (value < RF_CLASS_MIN_FIXED)
    {
        return 0U;
    }
    if (value >= RF_CLASS_MAX_FIXED)
    {
        return (RF_Class_t)(RF_NUM_CLASSES - 1U);
    }

    int32_t offset = (int32_t)value - (int32_t)RF_CLASS_MIN_FIXED;
    if (offset < 0)
    {
        offset = 0;
    }
    uint32_t cls = offset / (uint32_t)RF_CLASS_WIDTH_FIXED;

    if (cls >= RF_NUM_CLASSES)
    {
        cls = RF_NUM_CLASSES - 1U;
    }

    return (RF_Class_t)cls;
}

/**
 * @brief Sign of the difference of two class values (+1, -1 or 0).
 *
 * Used by RF_FlushResiduumRepeated() at the wrap of the repeated-
 * residue technique to decide whether and how a wrap point collapses.
 *
 * @param[in] hi First class value (minuend).
 * @param[in] lo Second class value (subtrahend).
 * @return +1 if hi>lo, -1 if hi<lo, 0 if hi==lo.
 */
static inline int32_t rf_sign_diff(RF_Class_t hi, RF_Class_t lo)
{
    const int32_t d = (int32_t)hi - (int32_t)lo;
    return (d > 0) ? 1 : ((d < 0) ? -1 : 0);
}

/**
 * @brief Push a classified turning point onto an arbitrary stack and
 *        check/close in cascade under the 4-point condition.
 *
 * Generic variant not bound to RF_Ctx_t::residuum, so the same
 * closure logic can be used for both the real residue stack and the
 * temporary stack of the repeated-residue technique
 * (see RF_FlushResiduumRepeated()).
 *
 * For each B-C closure the damage contribution from RF_DAMAGE_LUT
 * (indexed by class difference |B-C|) is added to *out_damage AND
 * the counter for that class difference in out_counts[|B-C|] is
 * incremented (range-pair histogram, see RF_GetRangePairCounts()).
 * B and C are removed from the stack by sliding D down two positions
 * onto B's slot.
 *
 * IMPORTANT PRECONDITION: the points passed in must (except the very
 * first) be true turning points — i.e. the sequence of push calls
 * must alternate (in the sense of rf_sign_diff()). If non-alternating
 * ("through") points are pushed, the 4-point condition can fire
 * falsely or miss a real closure (see RF_FlushResiduumRepeated()
 * history).
 *
 * Defensive plateau guard: if the new class_val equals the current
 * stack top exactly, it is NOT pushed (no-op, no closure needed — a
 * plateau carries no new direction information). Regularly this case
 * should never occur — the hysteresis filter categorically excludes
 * identical consecutive class values in the real residue, and
 * RF_FlushResiduumRepeated() 4-point-pushes running_extreme only when
 * it differs from residuum[count-1]. The check is therefore redundant,
 * but extra insurance against future callers or subtle call-order
 * bugs, without weakening the alternation precondition for true
 * turning points above.
 *
 * @param[in,out] stack      Stack array (caller-owned).
 * @param[in,out] count      Current fill of the stack.
 * @param[in]     capacity   Capacity of stack.
 * @param[in]     class_val  Class index of the new turning point.
 * @param[out]    out_damage Increased by all damage contributions
 *                            arising in this call (not overwritten).
 * @param[out]    out_counts Array of RF_NUM_CLASSES entries
 *                            (caller-owned) — out_counts[diff] is
 *                            incremented for each class difference
 *                            closed in this call (not overwritten).
 *                            May be NULL if the caller does not need
 *                            the range-pair histogram.
 *
 * @return RF_OK on success, RF_ERR_RESIDUUM_FULL if capacity is
 *         insufficient.
 */
static RF_Status_t rf_stack_push_and_close_buf(RF_Class_t *stack, uint16_t *count,
                                                uint16_t capacity, RF_Class_t class_val,
                                                uint32_t *out_damage, uint32_t *out_counts)
{
    if (*count >= capacity)
    {
        return RF_ERR_RESIDUUM_FULL;
    }

    if ((*count > 0U) && (stack[*count - 1U] == class_val))
    {
        /* Defensive plateau: same class value as the current stack
         * top — see function comment. Not a true turning point, so
         * no push and no closure check. */
        return RF_OK;
    }

    stack[(*count)++] = class_val;

    while (*count >= 4U)
    {
        const RF_Class_t a = stack[*count - 4U];
        const RF_Class_t b = stack[*count - 3U];
        const RF_Class_t c = stack[*count - 2U];
        const RF_Class_t d = stack[*count - 1U];

        const RF_Class_t min_bc = (b < c) ? b : c;
        const RF_Class_t max_bc = (b > c) ? b : c;
        const RF_Class_t min_ad = (a < d) ? a : d;
        const RF_Class_t max_ad = (a > d) ? a : d;

        if ((min_bc < min_ad) || (max_bc > max_ad))
        {
            break; /* B-C is not fully inside A-D: done */
        }

        const RF_Class_t diff = (b > c) ? (RF_Class_t)(b - c) : (RF_Class_t)(c - b);
        *out_damage += RF_DAMAGE_LUT[diff];
        if (out_counts != NULL)
        {
            out_counts[diff]++;
        }

        /* Remove B and C: D slides onto B's position, count -= 2 */
        stack[*count - 3U] = d;
        *count -= 2U;
    }

    return RF_OK;
}

/**
 * @brief Push a classified turning point onto the context residue
 *        stack and close in cascade (thin wrapper around
 *        rf_stack_push_and_close_buf()).
 *
 * @param[in,out] ctx        Rainflow context.
 * @param[in]     class_val  Class index of the new turning point.
 * @param[out]    out_damage Increased by all damage contributions
 *                            arising in this call (not overwritten).
 *
 * @return RF_OK on success, RF_ERR_RESIDUUM_FULL on a configuration
 *         error (see RF_MAX_RESIDUUM in rainflow.h) — unreachable in
 *         normal operation.
 */
static RF_Status_t rf_stack_push_and_close(RF_Ctx_t *ctx, RF_Class_t class_val,
                                            uint32_t *out_damage)
{
    /* Pass ctx->rp_counts through as out_counts: every closure from
     * RF_ProcessSample() is by definition an actual (not merely
     * predicted) count, same commit semantics as
     * ctx->accumulated_damage. */
    return rf_stack_push_and_close_buf(ctx->residuum, &ctx->residuum_count,
                                        RF_MAX_RESIDUUM, class_val, out_damage,
                                        ctx->rp_counts);
}

/**
 * @brief Add a single damage increment (max uint32_t) to a 96-bit
 *        accumulator (RF_Damage96_t), including carry from the low
 *        64-bit half (lo) into the high 32-bit half (hi).
 *
 * Overflow detection uses the standard unsigned-integer idiom: after
 * adding, the result is smaller than the original value only if a
 * carry occurred (unsigned arithmetic in C is modulo 2^N, so no
 * undefined behaviour). Overflow of hi itself (> 2^96 - 1 total) is
 * practically impossible for a single uint32_t increment per call and
 * is deliberately not handled here.
 *
 * @param[in,out] acc       96-bit accumulator.
 * @param[in]     increment Damage increment to add.
 */
static void rf_damage96_add(RF_Damage96_t *acc, uint32_t increment)
{
    const uint64_t old_lo = acc->lo;

    acc->lo += (uint64_t)increment;

    if (acc->lo < old_lo)
    {
        acc->hi++;
    }
}

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

RF_Status_t RF_Init(RF_Ctx_t *ctx, RF_Value_t hysteresis)
{
    if (ctx == NULL)
    {
        return RF_ERR_NULL_PTR;
    }

    /* ctx is deliberately zeroed even on error (rather than left
     * undefined) — a defined, if invalid, state is safer than
     * uninitialized memory if the return value is accidentally ignored. */
    memset(ctx, 0, sizeof(*ctx));
    ctx->hysteresis = hysteresis;
    ctx->slope = RF_SLOPE_UNKNOWN;
    /* ctx->accumulated_damage (lo and hi) is already 0 via memset(). */

    if (hysteresis < RF_CLASS_WIDTH_FIXED)
    {
        /* See the function comment in rainflow.h for the full
         * rationale: without this precondition the defensive plateau
         * guard in rf_stack_push_and_close_buf() can discard true
         * turning points and thus break the residue alternation
         * required for correctness. */
        return RF_ERR_INVALID_CONFIG;
    }

    return RF_OK;
}

RF_Status_t RF_ProcessSample(RF_Ctx_t *ctx, RF_Value_t sample,
                              uint32_t *out_damage_increment)
{
    if (ctx == NULL)
    {
        return RF_ERR_NULL_PTR;
    }

    /* A new sample reopens the stream: a previous repeated flush
     * must no longer count as already applied. */
    ctx->repeated_applied = false;

    /* Local working variable instead of writing *out_damage_increment
     * directly, because out_damage_increment is optional (NULL allowed)
     * — ctx->accumulated_damage must still be updated in every case. */
    uint32_t damage_increment = 0U;

    RF_Value_t confirmed_pt;
    const bool has_turning_point = rf_hysteresis_filter(ctx, sample, &confirmed_pt);

    if (has_turning_point)
    {
        const RF_Class_t class_val = rf_classify(confirmed_pt);
        const RF_Status_t st = rf_stack_push_and_close(ctx, class_val, &damage_increment);

        rf_damage96_add(&ctx->accumulated_damage, damage_increment);

        if (out_damage_increment != NULL)
        {
            *out_damage_increment = damage_increment;
        }

        return st;
    }

    if (out_damage_increment != NULL)
    {
        *out_damage_increment = 0U;
    }

    return RF_OK;
}

uint16_t RF_GetResiduumCount(const RF_Ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return 0U;
    }

    return ctx->residuum_count;
}

RF_Damage96_t RF_GetAccumulatedDamage(const RF_Ctx_t *ctx)
{
    if (ctx == NULL)
    {
        const RF_Damage96_t zero = { 0U, 0U };
        return zero;
    }

    return ctx->accumulated_damage;
}

double RF_GetAccumulatedDamageDouble(const RF_Ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return 0.0;
    }

    /* 2^64 — as a power of two exactly representable in double, so a
     * literal instead of math.h (pow()/ldexp()) — no extra libm
     * dependency. See the function comment in rainflow.h for the
     * precision note above 2^53. */
    return ((double)ctx->accumulated_damage.hi * 18446744073709551616.0)
           + (double)ctx->accumulated_damage.lo;
}

/**
 * @brief Close the currently open residue remainder with the ASTM
 *        E1049 repeated-residue technique: the residue is notionally
 *        appended to itself (as if this load segment repeated
 *        periodically) and the 4-point algorithm is run once over
 *        that. This finds extra, larger/nested full cycles that a
 *        simple pairwise count without a repetition assumption would
 *        miss. Works on a local temporary stack — ctx->residuum is
 *        unchanged during the calculation (see below for commit==true).
 *
 * Stage-1 candidate (running_extreme): the raw-signal candidate of
 * hysteresis stage 1 (see rf_hysteresis_filter()) is itself NOT yet
 * a confirmed turning point — a future sample could still extend it,
 * reverse it (inside the hysteresis band), or confirm it regularly.
 * For the cycle calculation here it is included anyway (as if the
 * signal had stopped right now): it represents the latest known
 * position of the raw signal, and a flush without it would
 * systematically ignore the last, possibly already more than one
 * class-width, branch. It is classified and — ONLY if it actually
 * represents a new point different from residuum[count-1] (i.e. a
 * direction has already been established, ctx->slope !=
 * RF_SLOPE_UNKNOWN, and its class differs from the last confirmed
 * point) — 4-point-pushed onto the stack (rfcnt feed_finalize /
 * interim). Only then is the sequence closure-free again; merely
 * appending without 4-point would leave inner cycles open whose
 * outer arm is the interim point, and would make the residue length
 * differ from rfcnt res_raw. Because the direction INTO
 * running_extreme is forced by the state machine to be opposite the
 * direction INTO residuum[count-1] (otherwise residuum[count-1]
 * would not yet have been confirmed), this point slots into the
 * alternation.
 *
 * Implementation: O(n) with constant extra work at the wrap (no
 * second pass, no buffer for an actually doubled sequence). Core
 * observation: after the interim 4-point, eff is a closure-free
 * stack AND internally guaranteed strictly alternating (the
 * hysteresis filter confirms a point only after a move > hysteresis
 * against the current direction — so a plateau can never occur
 * internally in the real residue). Therefore no push of the middle
 * points can trigger a closure (direct copy instead of individual
 * push calls), and a plateau can arise only at the single wrap where
 * the last point of the sequence meets the (notional) copy of its
 * first point.
 *
 * At that wrap three signs are compared (eff/eff_n denotes the
 * sequence possibly extended by running_extreme):
 *   d_prev = sign(eff[eff_n-1] - eff[eff_n-2])  (direction INTO eff[eff_n-1])
 *   d_wrap = sign(eff[0]       - eff[eff_n-1])  (direction ACROSS the wrap)
 *   d_next = sign(eff[1]       - eff[0])        (direction OUT OF eff[0])
 *
 * Case d_wrap != 0 (no value equality at the wrap):
 *   eff[eff_n-1] is dropped if d_wrap == d_prev (continuation of the
 *   same direction — eff[eff_n-1] is no longer a turning point).
 *   The copy of eff[0] is dropped if d_wrap == d_next (symmetric).
 *   These two cases are independent and can occur together.
 *
 * Case d_wrap == 0 (plateau at the wrap, eff[0] == eff[eff_n-1] —
 *   by the invariant above can occur ONLY here, never from real
 *   measured data): eff[eff_n-1] is always dropped (merged with the
 *   equal-valued copy of eff[0]). Whether the copy of eff[0] is also
 *   dropped depends on whether the merged point is itself still a
 *   true turning point:
 *     d_prev != d_next -> yes (simple merge, done).
 *     d_prev == d_next -> no (one step deeper, eff[1] becomes the
 *                         first point of the second copy).
 *   Proof that a third cascade step is NEVER needed: the next
 *   comparison direction would necessarily be -d_next (forced by the
 *   guaranteed internal alternation of the remaining sequence). A
 *   third collapse would need d_prev == -d_next; combined with the
 *   condition for the second step (d_prev == d_next) that yields
 *   d_prev == -d_prev, which is impossible for a sign other than 0.
 *   Verified exhaustively (n up to 9, 8 classes) and randomly (n up
 *   to 40, 64 classes) — see test_RF_FlushResiduumRepeated.c (there
 *   without the running_extreme extension; the extra sequence
 *   demonstrably satisfies the same alternation invariant and thus
 *   the same proof).
 *
 * The commit parameter decides whether the result is applied:
 *   commit == false (Predict): compute only, ctx remains completely
 *     unchanged in EVERY case (including running_extreme/slope).
 *   commit == true: the residue stays the 4-point-reduced sequence
 *     including the (possibly) adopted running_extreme — same as
 *     rfcnt res_raw after feed_finalize. Doubling only counts extra
 *     damage; it does not empty the stack. If running_extreme was
 *     included, it is treated as synthetically confirmed
 *     (ctx->slope = RF_SLOPE_UNKNOWN). A second flush without a new
 *     sample is a no-op (ctx->repeated_applied).
 *
 * @param[in,out] ctx                Rainflow context.
 * @param[in]     commit             true: update residue (and possibly
 *                                   stage-1 state, see above) after
 *                                   the calculation.
 *                                   false: compute only, leave ctx
 *                                   completely unchanged.
 * @param[out]    out_damage_increment Optional (may be NULL) —
 *                                     see rainflow.h. If given:
 *                                     summed damage from this
 *                                     calculation (not cumulative
 *                                     with previous calls).
 *
 * @return RF_OK on success, RF_ERR_NULL_PTR if ctx is NULL,
 *         RF_ERR_RESIDUUM_FULL on a configuration error
 *         (unreachable in normal operation).
 */
RF_Status_t RF_FlushResiduumRepeated(RF_Ctx_t *ctx, bool commit,
                                      uint32_t *out_damage_increment)
{
    if (ctx == NULL)
    {
        return RF_ERR_NULL_PTR;
    }

    if (commit && ctx->repeated_applied)
    {
        if (out_damage_increment != NULL)
        {
            *out_damage_increment = 0U;
        }
        return RF_OK;
    }
    /* Local working variable instead of writing *out_damage_increment
     * directly, because out_damage_increment is optional (NULL allowed)
     * — ctx->accumulated_damage (when commit == true) must still be
     * updated correctly in every case. */
    uint32_t damage_increment = 0U;

    /* Analogous local working array for range-pair counting — same
     * Predict safety as damage_increment: merged into ctx->rp_counts
     * only when commit == true. */
    uint32_t local_counts[RF_NUM_CLASSES];
    memset(local_counts, 0, sizeof(local_counts));

    /* Working copy: residue plus one slot for the interim. The
     * stage-1 candidate is 4-point-pushed onto the stack like the
     * rfcnt interim (not merely appended): otherwise inner cycles
     * stay open whose outer arm is the last candidate, and the
     * residue length differs from rfcnt res_raw. Before the push,
     * ctx->residuum may already hold RF_MAX_RESIDUUM confirmed
     * points; the extra slot takes the interim, then 4-point
     * reduces back to <= RF_MAX_RESIDUUM. */
    RF_Class_t eff[RF_MAX_RESIDUUM + 1U];
    uint16_t   eff_n = ctx->residuum_count;
    RF_Status_t st;

    if (eff_n > RF_MAX_RESIDUUM)
    {
        if (out_damage_increment != NULL)
        {
            *out_damage_increment = damage_increment;
        }
        return RF_ERR_RESIDUUM_FULL; /* configuration error, unreachable in normal operation */
    }
    memcpy(eff, ctx->residuum, (size_t)eff_n * sizeof(RF_Class_t));

    bool appended_running_extreme = false;

    if (ctx->has_running_extreme && (ctx->slope != RF_SLOPE_UNKNOWN))
    {
        const RF_Class_t re_class = rf_classify(ctx->running_extreme);

        if ((eff_n == 0U) || (re_class != eff[eff_n - 1U]))
        {
            st = rf_stack_push_and_close_buf(eff, &eff_n,
                    (uint16_t)(RF_MAX_RESIDUUM + 1U),
                    re_class, &damage_increment, local_counts);
            if (st != RF_OK)
            {
                if (out_damage_increment != NULL)
                {
                    *out_damage_increment = damage_increment;
                }
                return st;
            }
            appended_running_extreme = true;
        }
    }

    if (eff_n > RF_MAX_RESIDUUM)
    {
        /* 4-point did not bring the interim back to 2*RF_NUM_CLASSES
         * (including interim). */
        if (out_damage_increment != NULL)
        {
            *out_damage_increment = damage_increment;
        }
        return RF_ERR_RESIDUUM_FULL;
    }
    if (eff_n < 2U)
    {
        /* 0 or 1 point: nothing to double. Still carry over the
         * remainder if the interim point is the only new one. */
        if (commit)
        {
            memcpy(ctx->residuum, eff, (size_t)eff_n * sizeof(RF_Class_t));
            ctx->residuum_count = eff_n;
            rf_damage96_add(&ctx->accumulated_damage, damage_increment);
            for (size_t i = 0U; i < (size_t)RF_NUM_CLASSES; i++)
            {
                ctx->rp_counts[i] += local_counts[i];
            }
            if (appended_running_extreme)
            {
                ctx->slope = RF_SLOPE_UNKNOWN;
            }
            ctx->repeated_applied = true;
        }
        if (out_damage_increment != NULL)
        {
            *out_damage_increment = damage_increment;
        }
        return RF_OK;
    }

    const uint16_t n = eff_n;
    RF_Class_t local_stack[2U * (RF_MAX_RESIDUUM + 1U)];
    uint16_t local_count = n;

    /* Direct copy instead of a push loop: no push can trigger a
     * closure here, because after the interim 4-point, eff is already
     * a closure-free stack (see function comment). */
    memcpy(local_stack, eff, (size_t)n * sizeof(RF_Class_t));

    const int32_t d_prev = rf_sign_diff(eff[n - 1U], eff[n - 2U]);
    const int32_t d_wrap = rf_sign_diff(eff[0],      eff[n - 1U]);
    const int32_t d_next = rf_sign_diff(eff[1],      eff[0]);

    uint16_t skip_first;

    if (d_wrap == 0)
    {
        /* Wrap plateau (see function comment) — eff[n-1] is always
         * dropped; whether eff[0] is also dropped depends on the
         * cascade case proven to be at most one step. */
        local_count--;
        skip_first = (uint16_t)(d_prev == d_next);
    }
    else if (d_wrap == d_prev)
    {
        local_count--;
        skip_first = (uint16_t)(d_wrap == d_next);
    }
    else
    {
        skip_first = (uint16_t)(d_wrap == d_next);
    }

    for (uint16_t i = skip_first; i < n; i++)
    {
        st = rf_stack_push_and_close_buf(local_stack, &local_count,
                (uint16_t)(2U * (RF_MAX_RESIDUUM + 1U)), eff[i], &damage_increment,
                local_counts);

        if (st != RF_OK)
        {
            if (out_damage_increment != NULL)
            {
                *out_damage_increment = damage_increment;
            }
            return st; /* configuration error, unreachable in normal operation */
        }
    }

    /* Residuum = 4-point-reduced sequence including interim (rfcnt
     * res_raw). The doubling only counts additional damage.
     * If commit == false, ctx remains completely untouched. */
    if (commit)
    {
        memcpy(ctx->residuum, eff, (size_t)eff_n * sizeof(RF_Class_t));
        ctx->residuum_count = eff_n;
        rf_damage96_add(&ctx->accumulated_damage, damage_increment);

        for (size_t i = 0U; i < (size_t)RF_NUM_CLASSES; i++)
        {
            ctx->rp_counts[i] += local_counts[i];
        }

        if (appended_running_extreme)
        {
            /* running_extreme was treated synthetically as a confirmed
             * reversal point: reset the Level 1 state accordingly
             * as if it had just been confirmed,
             * with no known direction of movement—analogous to the
             * anchor case of the very first sample (see
             * rf_hysteresis_filter()). running_extreme itself remains
             * unchanged; its value already corresponds exactly to this
             * new anchor. */
            ctx->slope = RF_SLOPE_UNKNOWN;
        }
        ctx->repeated_applied = true;
    }

    if (out_damage_increment != NULL)
    {
        *out_damage_increment = damage_increment;
    }

    return RF_OK;
}

const uint32_t *RF_GetRangePairCounts(const RF_Ctx_t *ctx)
{
    if (ctx == NULL)
    {
        return NULL;
    }

    return ctx->rp_counts;
}
