/**
 * @file test_rainflow.c
 * @brief Standalone host tests for the rainflow engine (rainflow.h/.c).
 *
 * Pure assert()-based tests with no external framework, kept
 * deliberately simple (target system has no test infrastructure).
 * Compile and run, for example:
 *   gcc -std=c11 -Wall -Wextra -Wpedantic -I.. ../rainflow.c test_rainflow.c -o test_rainflow
 *   ./test_rainflow
 *
 * Supplies its own RF_DAMAGE_LUT (identity: LUT[i] == i) so expected
 * damage values can be checked by hand or script — rainflow_config.c
 * (the "real" LUT filled from the S-N curve) is deliberately NOT
 * linked in.
 *
 * Covers the public API against the CURRENT rainflow.h:
 *   - RF_Init() including RF_ERR_INVALID_CONFIG when
 *     hysteresis < RF_CLASS_WIDTH_FIXED
 *   - RF_ProcessSample() (optional out_damage_increment pointer)
 *   - RF_FlushResiduumRepeated() WITH commit parameter (Predict vs
 *     real flush) and inclusion of the still-unconfirmed stage-1
 *     candidate (running_extreme)
 *   - RF_GetAccumulatedDamage()/-Double() (96-bit accumulator)
 *   - RF_GetRangePairCounts()
 *
 * Detail correctness of RF_FlushResiduumRepeated() itself (wrap
 * cascade handling of the repeated-residue technique) is NOT verified
 * here, but exhaustively/randomized in
 * test_RF_FlushResiduumRepeated.c against an independent reference
 * implementation — this file checks the public API end-to-end on a
 * few examples recomputed by hand (or by a companion Python
 * reference script).
 */
#include "rainflow.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

const uint32_t RF_DAMAGE_LUT[RF_NUM_CLASSES] = {
    /* LUT[i] = i — makes expected test values checkable by hand. */
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
    96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
};

/* Q28.4 helper for readable test values in "real" units. */
static RF_Value_t q(float real_value)
{
    return (RF_Value_t)(real_value * (float)(1 << RF_FIXED_SHIFT));
}

/**
 * @brief Return the Q28.4 raw value that classifies EXACTLY to class k
 *        (k = 0 .. RF_NUM_CLASSES-1).
 *
 * rf_classify() rounds via (offset + RF_CLASS_WIDTH_FIXED/2) / RF_CLASS_WIDTH_FIXED;
 * for offset == k * RF_CLASS_WIDTH_FIXED integer division returns
 * exactly k (recomputed/verified in the companion Python reference
 * script). Lets test sequences be expressed directly in class indices,
 * independent of the actual numeric values of RF_CLASS_MIN_FIXED /
 * RF_CLASS_WIDTH_FIXED.
 */
static RF_Value_t class_center(int k)
{
    return (RF_Value_t)(RF_CLASS_MIN_FIXED + (int32_t)k * (int32_t)RF_CLASS_WIDTH_FIXED);
}

/* Comparison helper for the 96-bit damage accumulator against an
 * expected value that fits in a uint64_t (sufficient for all tests
 * here — a 96-bit overflow test would need billions of calls and is
 * not practical here). */
static bool damage96_equals_u64(RF_Damage96_t d, uint64_t expected)
{
    return (d.hi == 0U) && (d.lo == expected);
}

/* ------------------------------------------------------------------ */
/* Test: NULL-pointer handling and invalid configuration               */
/* ------------------------------------------------------------------ */
static void test_null_ptr_and_invalid_config(void)
{
    uint32_t damage = 0U;

    assert(RF_Init(NULL, RF_CLASS_WIDTH_FIXED) == RF_ERR_NULL_PTR);

    RF_Ctx_t ctx;
    assert(RF_Init(&ctx, RF_CLASS_WIDTH_FIXED) == RF_OK); /* boundary: exactly RF_CLASS_WIDTH_FIXED is still valid */

    /* Precondition hysteresis >= RF_CLASS_WIDTH_FIXED is strictly enforced. */
    RF_Ctx_t ctx_invalid;
    assert(RF_Init(&ctx_invalid, (RF_Value_t)(RF_CLASS_WIDTH_FIXED - 1)) == RF_ERR_INVALID_CONFIG);

    assert(RF_ProcessSample(NULL, q(1.0f), &damage) == RF_ERR_NULL_PTR);
    /* out_damage_increment is optional (unlike an earlier API version)
     * — NULL is not an error here as long as ctx is valid. */
    assert(RF_ProcessSample(&ctx, q(1.0f), NULL) == RF_OK);

    assert(RF_FlushResiduumRepeated(NULL, true, &damage) == RF_ERR_NULL_PTR);
    assert(RF_FlushResiduumRepeated(&ctx, true, NULL) == RF_OK); /* also optional */

    assert(RF_GetResiduumCount(NULL) == 0U);
    assert(damage96_equals_u64(RF_GetAccumulatedDamage(NULL), 0U));
    assert(RF_GetAccumulatedDamageDouble(NULL) == 0.0);
    assert(RF_GetRangePairCounts(NULL) == NULL);

    printf("PASS: test_null_ptr_and_invalid_config\n");
}

/* ------------------------------------------------------------------ */
/* Test: stage 1 — first raw sample is always confirmed as the anchor  */
/* ------------------------------------------------------------------ */
static void test_first_sample_is_anchor(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);

    uint32_t damage = 0U;
    RF_Status_t st = RF_ProcessSample(&ctx, class_center(70), &damage);

    assert(st == RF_OK);
    assert(damage == 0U); /* count < 4, no closure possible */
    assert(RF_GetResiduumCount(&ctx) == 1U);
    assert(damage96_equals_u64(RF_GetAccumulatedDamage(&ctx), 0U));

    printf("PASS: test_first_sample_is_anchor\n");
}

/* ------------------------------------------------------------------ */
/* Test: small reversals inside the hysteresis band are ignored        */
/* ------------------------------------------------------------------ */
static void test_hysteresis_suppresses_noise(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED); /* hysteresis = RF_CLASS_WIDTH_FIXED (minimum) */

    uint32_t damage = 0U;
    RF_ProcessSample(&ctx, class_center(60), &damage);                     /* anchor */
    RF_ProcessSample(&ctx, (RF_Value_t)(class_center(60) + 2000), &damage); /* RISING direction established */
    uint32_t inc = 0U;
    RF_ProcessSample(&ctx, (RF_Value_t)(class_center(60) + 2000 - 100), &inc); /* reversal 100 < 500 hysteresis: ignored */

    assert(inc == 0U);
    assert(RF_GetResiduumCount(&ctx) == 1U); /* only the anchor is in the residue so far */

    printf("PASS: test_hysteresis_suppresses_noise\n");
}

/* ------------------------------------------------------------------ */
/* Test: a monotonically rising signal produces no cycles              */
/* ------------------------------------------------------------------ */
static void test_monotonic_signal_produces_no_cycles(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);

    uint32_t total_damage = 0U;
    for (int k = 0; k <= 10; k++)
    {
        uint32_t inc = 0U;
        RF_ProcessSample(&ctx, class_center(k * 10), &inc);
        total_damage += inc;
    }

    assert(total_damage == 0U);
    assert(RF_GetResiduumCount(&ctx) == 1U); /* only the start anchor */

    printf("PASS: test_monotonic_signal_produces_no_cycles\n");
}

/* ------------------------------------------------------------------ */
/* Test: values below RF_CLASS_MIN_FIXED are clamped to class 0        */
/*                                                                      */
/* Note: the mirrored test for the UPPER bound is not meaningful with  */
/* the current rainflow_config.h (placeholder values): RF_CLASS_MAX    */
/* sits at 2500.0 real units, while the largest representable Q28.4    */
/* raw value (RF_Value_t, int32_t) is only about 2047.94 — the upper   */
/* clamp clause in rf_classify() is therefore unreachable by any valid */
/* RF_Value_t with this configuration (a configuration question, see   */
/* the TODO in rainflow_config.h — not a bug in rf_classify() itself). */
/* ------------------------------------------------------------------ */
static void test_out_of_range_values_are_clamped_to_lowest_class(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);

    uint32_t damage = 0U;
    RF_ProcessSample(&ctx, class_center(60), &damage); /* anchor, class 60 */
    RF_ProcessSample(&ctx, q(-2000.0f), &damage);       /* < RF_CLASS_MIN_FIXED (-1500): only FALLING direction established */
    assert(RF_GetResiduumCount(&ctx) == 1U);

    uint32_t inc = 0U;
    RF_ProcessSample(&ctx, class_center(60), &inc); /* reversal: confirms the -2000 point as a valley */
    assert(RF_GetResiduumCount(&ctx) == 2U);
    assert(ctx.residuum[1] == 0U); /* clamped to lowest class */

    printf("PASS: test_out_of_range_values_are_clamped_to_lowest_class\n");
}

/* ------------------------------------------------------------------ */
/* Test: full pipeline including cascade during the stream, Predict    */
/* flush (commit=false) and subsequent real flush (commit=true) —      */
/* including inclusion of the still-unconfirmed stage-1 candidate      */
/* (running_extreme).                                                  */
/*                                                                      */
/* Input sequence as class indices (driven directly via                */
/* class_center()): 50 (anchor), 60, 55, 70, 45, 44, 58, 58 (repeat).  */
/* Expected values were recomputed with an independent Python replica  */
/* of rf_classify()/rf_stack_push_and_close_buf()/                     */
/* RF_FlushResiduumRepeated() (identical logic, separately             */
/* implemented):                                                       */
/*   - during the stream exactly one cascade closure: when confirming  */
/*     class 70, (60,55) closes, diff=5                                */
/*   - afterwards residue = [50, 70, 44], running_extreme candidate    */
/*     = 58 (RISING direction, still NOT confirmed)                    */
/*   - Predict flush (commit=false): damage increment 34, ctx remains  */
/*     completely unchanged                                            */
/*   - real flush (commit=true): the same increment 34 is added to     */
/*     ctx->accumulated_damage; residue collapses to the last point    */
/*     [58], running_extreme candidate is treated as synthetically     */
/*     confirmed -> stage-1 state reset                                */
/*   - range-pair histogram after the real flush: exactly one entry    */
/*     each at class differences 5, 8 and 26                           */
/*   - total damage over the whole test: 5 + 34 = 39                   */
/* ------------------------------------------------------------------ */
static void test_full_pipeline_cascade_predict_and_commit_flush(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);

    const int ks[] = { 50, 60, 55, 70, 45, 44, 58, 58 };
    const size_t n = sizeof(ks) / sizeof(ks[0]);
    uint32_t total_damage = 0U;

    for (size_t i = 0U; i < n; i++)
    {
        uint32_t inc = 0U;
        RF_Status_t st = RF_ProcessSample(&ctx, class_center(ks[i]), &inc);
        assert(st == RF_OK);
        total_damage += inc;
    }

    assert(total_damage == 5U); /* the one cascade closure during the stream (60,55) */
    assert(RF_GetResiduumCount(&ctx) == 3U);
    assert(ctx.residuum[0] == 50U);
    assert(ctx.residuum[1] == 70U);
    assert(ctx.residuum[2] == 44U);
    assert(damage96_equals_u64(RF_GetAccumulatedDamage(&ctx), 5U));

    /* Predict flush: computes the increment including running_extreme
     * (58, RISING, still unconfirmed), but does NOT change ctx. */
    RF_Ctx_t ctx_snapshot_before_predict = ctx;
    uint32_t predict_damage = 0U;
    assert(RF_FlushResiduumRepeated(&ctx, false, &predict_damage) == RF_OK);
    assert(predict_damage == 34U);
    assert(memcmp(&ctx, &ctx_snapshot_before_predict, sizeof(ctx)) == 0); /* truly unchanged */

    /* Real flush: same increment, this time applied. */
    uint32_t commit_damage = 0U;
    assert(RF_FlushResiduumRepeated(&ctx, true, &commit_damage) == RF_OK);
    assert(commit_damage == 34U);
    total_damage += commit_damage;

    assert(RF_GetResiduumCount(&ctx) == 1U);
    assert(ctx.residuum[0] == 58U); /* last point of the extended sequence remains as the base */
    assert(total_damage == 39U);
    assert(damage96_equals_u64(RF_GetAccumulatedDamage(&ctx), 39U));
    assert(RF_GetAccumulatedDamageDouble(&ctx) == 39.0);

    const uint32_t *rp = RF_GetRangePairCounts(&ctx);
    assert(rp != NULL);
    for (size_t diff = 0U; diff < RF_NUM_CLASSES; diff++)
    {
        const uint32_t expected = ((diff == 5U) || (diff == 8U) || (diff == 26U)) ? 1U : 0U;
        assert(rp[diff] == expected);
    }

    printf("PASS: test_full_pipeline_cascade_predict_and_commit_flush\n");
}

/* ------------------------------------------------------------------ */
/* Test: after a real flush with no remaining stage-1 candidate        */
/* (slope == RF_SLOPE_UNKNOWN, see previous test) a further flush      */
/* call yields 0 — there is nothing left to close (only a single       */
/* point in the residue, no running_extreme that can be included).     */
/* ------------------------------------------------------------------ */
static void test_flush_after_full_commit_yields_zero(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);

    const int ks[] = { 50, 60, 55, 70, 45, 44, 58, 58 };
    for (size_t i = 0U; i < sizeof(ks) / sizeof(ks[0]); i++)
    {
        RF_ProcessSample(&ctx, class_center(ks[i]), NULL);
    }
    RF_FlushResiduumRepeated(&ctx, true, NULL);

    RF_Ctx_t ctx_before = ctx;
    uint32_t damage = 123U; /* canary value — must be overwritten to 0 */
    assert(RF_FlushResiduumRepeated(&ctx, true, &damage) == RF_OK);
    assert(damage == 0U);
    assert(memcmp(&ctx, &ctx_before, sizeof(ctx)) == 0); /* nothing to do -> nothing changed */

    printf("PASS: test_flush_after_full_commit_yields_zero\n");
}

/* ------------------------------------------------------------------ */
/* Test: flush right after the very first sample (only the anchor in   */
/* the residue, no stage-1 candidate with a known direction) also      */
/* yields 0 — a single point alone can never form a cycle.             */
/* ------------------------------------------------------------------ */
static void test_predict_flush_with_single_point_yields_zero(void)
{
    RF_Ctx_t ctx;
    RF_Init(&ctx, RF_CLASS_WIDTH_FIXED);
    RF_ProcessSample(&ctx, class_center(70), NULL);

    uint32_t damage = 0U;
    assert(RF_FlushResiduumRepeated(&ctx, false, &damage) == RF_OK);
    assert(damage == 0U);
    assert(RF_GetResiduumCount(&ctx) == 1U); /* Predict: unchanged */

    printf("PASS: test_predict_flush_with_single_point_yields_zero\n");
}

int main(void)
{
    test_null_ptr_and_invalid_config();
    test_first_sample_is_anchor();
    test_hysteresis_suppresses_noise();
    test_monotonic_signal_produces_no_cycles();
    test_out_of_range_values_are_clamped_to_lowest_class();
    test_full_pipeline_cascade_predict_and_commit_flush();
    test_flush_after_full_commit_yields_zero();
    test_predict_flush_with_single_point_yields_zero();

    printf("\nALL TESTS PASSED\n");
    return 0;
}
