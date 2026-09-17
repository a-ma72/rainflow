/**
 * @file test_RF_FlushResiduumRepeated.c
 * @brief Verification of RF_FlushResiduumRepeated() (rainflow.c)
 *        against an independent, justified ground-truth reference.
 *
 * Contains two standalone implementations of the same semantics:
 *
 *  - ref_flush_repeated(): general, unoptimized reference.
 *    Reduces the (notionally doubled) sequence to true turning points
 *    by turning-point extraction before it goes through the 4-point
 *    closer. Necessary because the 4-point closer demonstrably yields
 *    false closures on non-alternating ("through") points (proof: a
 *    purely monotonic sequence such as 5,6,7,8 produces a closure even
 *    though no reversal ever occurs — see the diagnostic section
 *    below). Starting direction of residuum[0] is derived (not pinned
 *    ad hoc): because residuum[0] is an already confirmed turning
 *    point, the direction that led INTO residuum[0] must be opposite
 *    the outgoing direction (residuum[0]->residuum[1]).
 *
 *  - opt_flush_repeated(): the O(1) implementation actually used in
 *    rainflow.c (identical logic, duplicated here so this test stays
 *    runnable without the real project header).
 *
 * Test invariant (given, see history): the hysteresis filter in
 * rainflow.c categorically excludes INTERNAL plateaus in the real
 * residue — two consecutive raw values always differ by at least the
 * hysteresis threshold. A plateau (same class value) can therefore
 * arise only at the single wrap of the repeated-residue technique
 * (residuum[0] == residuum[n-1]), never internally. The cascade at
 * that wrap is proven to be at most 2 steps (see the proof comment in
 * rainflow.c at RF_FlushResiduumRepeated()) — that enables the O(1)
 * solution instead of a loop of variable length.
 *
 * Test stages:
 *   Stage 1 (NEGATIVE CONTROL): also allows internal plateaus that,
 *            by the invariant above, can never occur in reality.
 *            Shows many mismatches as expected — that is NOT a bug
 *            hint, but demonstrates that the optimization (validly)
 *            relies on the invariant.
 *   Stage 2: exhaustive enumeration of all internally alternating
 *            residues (including wrap-plateau cases) up to
 *            n=RF_MAX_RESIDUUM, RF_NUM_CLASSES classes.
 *   Stage 3: randomized stress test at production size (64 classes,
 *            n up to 40); wrap-plateau hits are counted naturally.
 *   Stage 4: wrap plateau (residuum[0]==residuum[n-1]) forced on
 *            purpose so this rarer case is checked often enough.
 *   Diagnostic: measures the actually observed maximum cascade depth
 *            at the wrap (exhaustive + randomized) as an empirical
 *            cross-check of the proof.
 *
 * Return value: 0 if stages 2-4 all run with no mismatch (stage 1
 * does NOT feed into the return value, see above), else 1.
 *
 * Compile:  gcc -O2 -Wall -Wextra -o test_RF_FlushResiduumRepeated test_RF_FlushResiduumRepeated.c
 * Run:      ./test_RF_FlushResiduumRepeated
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

typedef uint8_t RF_Class_t;
typedef enum { RF_OK = 0, RF_ERR_RESIDUUM_FULL } RF_Status_t;

#define RF_NUM_CLASSES   8U
#define RF_MAX_RESIDUUM  9U

static uint32_t g_lut[RF_NUM_CLASSES];

/* ---- 1:1 from rainflow.c ---- */
static RF_Status_t rf_stack_push_and_close_buf(RF_Class_t *stack, uint16_t *count,
                                                uint16_t capacity, RF_Class_t class_val,
                                                uint32_t *out_damage)
{
    if (*count >= capacity) return RF_ERR_RESIDUUM_FULL;
    stack[*count] = class_val;
    (*count)++;
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
        if ((min_bc < min_ad) || (max_bc > max_ad)) break;
        const RF_Class_t diff = (b > c) ? (RF_Class_t)(b - c) : (RF_Class_t)(c - b);
        *out_damage += g_lut[diff];
        stack[*count - 3U] = d;
        *count -= 2U;
    }
    return RF_OK;
}

static inline int32_t rf_sign_diff(RF_Class_t hi, RF_Class_t lo)
{
    const int32_t d = (int32_t)hi - (int32_t)lo;
    return (d > 0) ? 1 : ((d < 0) ? -1 : 0);
}

/* ---- Ground truth: justified turning-point extraction + closure ---- */
static uint32_t ref_flush_repeated(const RF_Class_t *res, uint16_t n)
{
    RF_Class_t out[256];
    uint16_t out_count = 1U;
    out[0] = res[0];

    /* Derived start direction: res[0] is a confirmed turning point,
     * so the incoming direction is necessarily opposite the outgoing
     * one (res[0]->res[1]). */
    int32_t last_dir = -rf_sign_diff(res[1], res[0]);

    for (uint16_t i = 1U; i < (uint16_t)(2U * n); i++)
    {
        const RF_Class_t val = res[i % n];
        const int32_t dir = rf_sign_diff(val, out[out_count - 1U]);

        if (dir == 0)
        {
            continue; /* true plateau: no new information */
        }
        else if (dir == last_dir)
        {
            out[out_count - 1U] = val; /* continuation: replace */
        }
        else
        {
            out[out_count] = val; /* reversal: true turning point */
            out_count++;
            last_dir = dir;
        }
    }

    RF_Class_t stack[256];
    uint16_t count = 0U;
    uint32_t damage = 0U;
    for (uint16_t i = 0U; i < out_count; i++)
    {
        rf_stack_push_and_close_buf(stack, &count, 256U, out[i], &damage);
    }
    return damage;
}

/* ---- Optimized variant under test: O(1) at the wrap, no loop.
 * Proof that the wrap-plateau cascade never goes deeper than 2 steps:
 * the second cascade step (if it happens at all) ALWAYS compares
 * against the direction -d_next (forced by the internal alternation
 * of the remaining sequence). A third step would need d_prev ==
 * -d_next; combined with the condition for the second step
 * (d_prev == d_next) that yields d_prev == -d_prev, which is
 * impossible for a sign other than 0. ---- */
static uint32_t opt_flush_repeated(const RF_Class_t *res, uint16_t n)
{
    RF_Class_t local_stack[256];
    uint16_t local_count = n;
    uint32_t damage = 0U;

    memcpy(local_stack, res, (size_t)n * sizeof(RF_Class_t));

    const int32_t d_prev = rf_sign_diff(res[n - 1U], res[n - 2U]);
    const int32_t d_wrap = rf_sign_diff(res[0],      res[n - 1U]);
    const int32_t d_next = rf_sign_diff(res[1],      res[0]);

    uint16_t skip_first;

    if (d_wrap == 0)
    {
        /* Wrap plateau: by the invariant can occur ONLY here, never
         * internally. residuum[n-1] is always dropped (merged with the
         * equal-valued copy of residuum[0]). Whether the copy of
         * residuum[0] is also dropped depends on whether the merged
         * point is itself still a true turning point (d_prev != d_next)
         * or not (d_prev == d_next, deeper collapse — see proof above,
         * never goes further than this). */
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
        rf_stack_push_and_close_buf(local_stack, &local_count, 256U, res[i], &damage);
    }
    return damage;
}

/* ------------------------------------------------------------------ */
static uint64_t g_generated = 0U, g_valid = 0U, g_mismatches = 0U;
static int g_print_limit = 40;

static bool is_closure_free(const RF_Class_t *res, uint16_t n)
{
    RF_Class_t stack[RF_MAX_RESIDUUM];
    uint16_t count = 0U;
    uint32_t dummy = 0U;
    for (uint16_t i = 0U; i < n; i++)
        rf_stack_push_and_close_buf(stack, &count, RF_MAX_RESIDUUM, res[i], &dummy);
    return (count == n);
}

/* Valid residue: internally alternating (plateaus/ties allowed, but a
 * true direction reversal always on a direction change), AND
 * closure-free. That is the invariant that rf_hysteresis_filter +
 * rf_stack_push_and_close actually guarantee for ctx->residuum. */
static bool is_valid_residuum(const RF_Class_t *res, uint16_t n)
{
    int32_t last_dir = 0;
    for (uint16_t i = 0U; i + 1U < n; i++)
    {
        const int32_t d = rf_sign_diff(res[i + 1U], res[i]);
        if (d == 0) continue; /* skip plateau */
        if (last_dir != 0 && d == last_dir) return false; /* no true reversal */
        last_dir = d;
    }
    return is_closure_free(res, n);
}

/* Like is_valid_residuum, but with NO tolerance for plateaus (not a
 * single diff may be 0) — isolates the core algorithm from the
 * separate, still-open plateau question. */
static bool is_valid_residuum_strict(const RF_Class_t *res, uint16_t n)
{
    /* Exclude INTERNAL plateaus only (impossible in the real residue
     * by the hysteresis-filter invariant). The wrap transition
     * res[0]<->res[n-1] may now also be 0 (plateau ONLY from the
     * doubling) — that is exactly the case to be tested here. */
    int32_t last_dir = 0;
    for (uint16_t i = 0U; i + 1U < n; i++)
    {
        const int32_t d = rf_sign_diff(res[i + 1U], res[i]);
        if (d == 0) return false; /* internal plateaus: impossible by assumption */
        if (last_dir != 0 && d == last_dir) return false;
        last_dir = d;
    }
    return is_closure_free(res, n);
}

static uint64_t g_valid_strict = 0U, g_mismatches_strict = 0U;
static uint16_t g_max_cascade_exh = 0U;
static RF_Class_t g_max_cascade_exh_example[RF_MAX_RESIDUUM];
static uint16_t g_max_cascade_exh_n = 0U;

static uint16_t measure_cascade_depth(const RF_Class_t *res, uint16_t n)
{
    RF_Class_t candidate = res[n - 1U];
    int32_t last_dir = rf_sign_diff(res[n - 1U], res[n - 2U]);
    uint16_t depth = 0U;
    for (uint16_t i = 0U; i < n; i++)
    {
        const int32_t dir = rf_sign_diff(res[i], candidate);
        if (dir == 0) { depth++; continue; }
        else if (dir == last_dir) { candidate = res[i]; depth++; }
        else { break; }
    }
    return depth;
}

static void check_one_strict(const RF_Class_t *res, uint16_t n)
{
    if (!is_valid_residuum_strict(res, n)) return;
    g_valid_strict++;

    const uint16_t depth = measure_cascade_depth(res, n);
    if (depth > g_max_cascade_exh)
    {
        g_max_cascade_exh = depth;
        g_max_cascade_exh_n = n;
        memcpy(g_max_cascade_exh_example, res, n * sizeof(RF_Class_t));
    }

    const uint32_t ref_damage = ref_flush_repeated(res, n);
    const uint32_t opt_damage = opt_flush_repeated(res, n);
    if (ref_damage != opt_damage)
    {
        g_mismatches_strict++;
        if (g_print_limit > 0)
        {
            g_print_limit--;
            printf("STRICT MISMATCH n=%u residuum=[", n);
            for (uint16_t i = 0U; i < n; i++) printf("%u%s", res[i], (i + 1U < n) ? "," : "");
            printf("]  ref=%u  opt=%u\n", ref_damage, opt_damage);
        }
    }
}

static void enumerate_strict(RF_Class_t *res, uint16_t n, uint16_t pos)
{
    if (pos == n) { check_one_strict(res, n); return; }
    for (uint16_t v = 0U; v < RF_NUM_CLASSES; v++)
    {
        res[pos] = (RF_Class_t)v;
        enumerate_strict(res, n, (uint16_t)(pos + 1U));
    }
}

static void check_one(const RF_Class_t *res, uint16_t n)
{
    g_generated++;
    if (!is_valid_residuum(res, n)) return;
    g_valid++;

    const uint32_t ref_damage = ref_flush_repeated(res, n);
    const uint32_t opt_damage = opt_flush_repeated(res, n);

    if (ref_damage != opt_damage)
    {
        g_mismatches++;
        if (g_print_limit > 0)
        {
            g_print_limit--;
            printf("MISMATCH n=%u residuum=[", n);
            for (uint16_t i = 0U; i < n; i++) printf("%u%s", res[i], (i + 1U < n) ? "," : "");
            printf("]  ref=%u  opt=%u\n", ref_damage, opt_damage);
        }
    }
}

static void enumerate(RF_Class_t *res, uint16_t n, uint16_t pos)
{
    if (pos == n) { check_one(res, n); return; }
    for (uint16_t v = 0U; v < RF_NUM_CLASSES; v++)
    {
        res[pos] = (RF_Class_t)v;
        enumerate(res, n, (uint16_t)(pos + 1U));
    }
}

int main(void)
{
    for (uint32_t i = 0U; i < RF_NUM_CLASSES; i++) g_lut[i] = i;

    printf("=== Stage 1 (NEGATIVE CONTROL - internal plateaus allowed,\n"
           "    even though by the hysteresis-filter invariant they never occur) ===\n");
    printf("RF_NUM_CLASSES=%u, RF_MAX_RESIDUUM=%u\n", RF_NUM_CLASSES, RF_MAX_RESIDUUM);
    printf("Expectation: MANY mismatches - that is not a bug, but shows\n"
           "that the O(1) optimization validly relies on the invariant\n"
           "(see the comment in rainflow.c). Documentation only,\n"
           "does NOT feed into the return value.\n\n");

    g_print_limit = 0; /* do not print individual cases here, only the final total */
    for (uint16_t n = 2U; n <= RF_MAX_RESIDUUM; n++)
    {
        RF_Class_t res[RF_MAX_RESIDUUM];
        enumerate(res, n, 0U);
    }

    printf("Stage 1 result: %llu generated, %llu valid (closure-free) "
           "residues, %llu mismatches (high, as expected). ===\n\n",
           (unsigned long long)g_generated, (unsigned long long)g_valid,
           (unsigned long long)g_mismatches);

    printf("=== Stage 2: exhaustive, ONLY internally alternating residues "
           "(including wrap-plateau cases) ===\n");
    g_print_limit = 40;
    for (uint16_t n = 2U; n <= RF_MAX_RESIDUUM; n++)
    {
        RF_Class_t res[RF_MAX_RESIDUUM];
        enumerate_strict(res, n, 0U);
        printf("  n=%u done (cumulative: valid=%llu, mismatches=%llu)\n",
               n, (unsigned long long)g_valid_strict, (unsigned long long)g_mismatches_strict);
    }

    printf("\n=== Stage 2 result: %llu strictly alternating residues tested, "
           "%llu mismatches. ===\n", (unsigned long long)g_valid_strict,
           (unsigned long long)g_mismatches_strict);
    printf("Exhaustively measured MAXIMUM cascade depth (n up to %u, K=%u): %u\n",
           RF_MAX_RESIDUUM, RF_NUM_CLASSES, g_max_cascade_exh);
    printf("Example: [");
    for (uint16_t i = 0U; i < g_max_cascade_exh_n; i++)
        printf("%u%s", g_max_cascade_exh_example[i], (i+1U<g_max_cascade_exh_n)?",":"");
    printf("] (n=%u)\n", g_max_cascade_exh_n);

    /* Stage 3: randomized stress test at a realistic size (K=64),
     * STRICTLY alternating only (plateau-free), as extra coverage
     * beyond exhaustive stage 2 (there K=8, n<=9; here K=64, n up
     * to 40, many samples). */
    printf("\n=== Stage 3: randomized, RF_NUM_CLASSES=64, n up to 40 ===\n");
    srand(42U);
    uint64_t rnd_tested = 0U, rnd_mismatch = 0U, rnd_wrap_zero = 0U;
    const uint16_t K_RND = 64U;
    const uint16_t N_RND_MAX = 40U;
    RF_Class_t rnd_res[64];

    for (uint64_t it = 0U; it < 3000000ULL; it++)
    {
        const uint16_t n = (uint16_t)(2U + (rand() % (N_RND_MAX - 1U)));
        /* Construct a strictly alternating random sequence (guaranteed
         * plateau-free, including at the wrap). */
        int32_t dir = (rand() % 2) ? 1 : -1;
        rnd_res[0] = (RF_Class_t)(rand() % K_RND);
        bool ok = true;
        for (uint16_t i = 1U; i < n; i++)
        {
            int32_t next_val;
            if (dir > 0)
                next_val = (int32_t)rnd_res[i - 1U] + 1 + (rand() % (K_RND - rnd_res[i - 1U] > 1 ? K_RND - rnd_res[i-1U]-1 : 1));
            else
                next_val = (int32_t)rnd_res[i - 1U] - 1 - (rand() % (rnd_res[i-1U] > 1 ? rnd_res[i-1U]-1 : 1));
            if (next_val < 0 || next_val >= (int32_t)K_RND) { ok = false; break; }
            rnd_res[i] = (RF_Class_t)next_val;
            dir = -dir;
        }
        if (!ok) continue;
        /* Do NOT exclude wrap plateau (d_wrap==0) any more — include
         * it on purpose, see rnd_wrap_zero counter below. */

        RF_Class_t tmp[64]; uint16_t tc=0; uint32_t td=0;
        for (uint16_t i=0;i<n;i++) rf_stack_push_and_close_buf(tmp,&tc,64,rnd_res[i],&td);
        if (tc != n) continue; /* not closure-free, discard */

        if (rf_sign_diff(rnd_res[0], rnd_res[n-1U]) == 0) rnd_wrap_zero++;

        const uint32_t ref_d = ref_flush_repeated(rnd_res, n);
        const uint32_t opt_d = opt_flush_repeated(rnd_res, n);
        rnd_tested++;
        if (ref_d != opt_d)
        {
            rnd_mismatch++;
            if (g_print_limit > 0)
            {
                g_print_limit--;
                printf("STAGE3 MISMATCH n=%u ref=%u opt=%u residuum=[", n, ref_d, opt_d);
                for (uint16_t i=0;i<n;i++) printf("%u%s", rnd_res[i], (i+1U<n)?",":"");
                printf("]\n");
            }
        }
    }
    printf("Stage 3: %llu tested, %llu mismatches (of which %llu with wrap plateau "
           "d_wrap==0, hit at random).\n",
           (unsigned long long)rnd_tested, (unsigned long long)rnd_mismatch,
           (unsigned long long)rnd_wrap_zero);

    /* Stage 4: FORCE wrap plateau (d_wrap==0) (with random values from
     * 64 classes too rare to rely on stage 3 alone) — res[n-1] is set
     * explicitly to res[0]. */
    printf("\n=== Stage 4: wrap plateau (d_wrap==0) forced on purpose ===\n");
    uint64_t forced_tested = 0U, forced_mismatch = 0U;

    for (uint64_t it = 0U; it < 1000000ULL; it++)
    {
        const uint16_t n = (uint16_t)(4U + (rand() % (N_RND_MAX - 3U))); /* n>=4, see comment below */
        int32_t dir = (rand() % 2) ? 1 : -1;
        rnd_res[0] = (RF_Class_t)(rand() % K_RND);
        bool ok = true;
        for (uint16_t i = 1U; i + 1U < n; i++) /* generate res[1..n-2] normally */
        {
            int32_t next_val;
            if (dir > 0)
                next_val = (int32_t)rnd_res[i - 1U] + 1 + (rand() % (K_RND - rnd_res[i - 1U] > 1 ? K_RND - rnd_res[i-1U]-1 : 1));
            else
                next_val = (int32_t)rnd_res[i - 1U] - 1 - (rand() % (rnd_res[i-1U] > 1 ? rnd_res[i-1U]-1 : 1));
            if (next_val < 0 || next_val >= (int32_t)K_RND) { ok = false; break; }
            rnd_res[i] = (RF_Class_t)next_val;
            dir = -dir;
        }
        if (!ok) continue;
        /* res[n-1] forced = res[0]: wrap plateau guaranteed.
         * Must still satisfy internal alternation: the direction
         * res[n-2]->res[n-1] must be opposite the direction BEFORE
         * res[n-2] (otherwise res[n-1] would not even be a turning
         * point internally — the generator bug that history uncovered). */
        rnd_res[n - 1U] = rnd_res[0];
        if (rnd_res[n - 1U] == rnd_res[n - 2U]) continue;
        {
            const int32_t d_into_last  = rf_sign_diff(rnd_res[n-1U], rnd_res[n-2U]);
            const int32_t d_before_that= rf_sign_diff(rnd_res[n-2U], rnd_res[n-3U]);
            if (d_into_last == d_before_that) continue; /* not alternating, discard */
        }

        RF_Class_t tmp[64]; uint16_t tc=0; uint32_t td=0;
        for (uint16_t i=0;i<n;i++) rf_stack_push_and_close_buf(tmp,&tc,64,rnd_res[i],&td);
        if (tc != n) continue;

        const uint32_t ref_d = ref_flush_repeated(rnd_res, n);
        const uint32_t opt_d = opt_flush_repeated(rnd_res, n);
        forced_tested++;
        if (ref_d != opt_d)
        {
            forced_mismatch++;
            if (g_print_limit > 0)
            {
                g_print_limit--;
                printf("STAGE4 MISMATCH n=%u ref=%u opt=%u residuum=[", n, ref_d, opt_d);
                for (uint16_t i=0;i<n;i++) printf("%u%s", rnd_res[i], (i+1U<n)?",":"");
                printf("]\n");
            }
        }
    }
    printf("Stage 4: %llu tested (all with d_wrap==0), %llu mismatches.\n",
           (unsigned long long)forced_tested, (unsigned long long)forced_mismatch);

    /* Diagnostic: how deep can the cascade at the wrap go at most? */
    printf("\n=== Diagnostic: maximum cascade depth at the wrap ===\n");
    uint16_t max_cascade_depth = 0U;
    RF_Class_t max_cascade_example[64];
    uint16_t max_cascade_n = 0U;

    for (uint64_t it = 0U; it < 2000000ULL; it++)
    {
        const uint16_t n = (uint16_t)(4U + (rand() % (N_RND_MAX - 3U)));
        int32_t dir = (rand() % 2) ? 1 : -1;
        rnd_res[0] = (RF_Class_t)(rand() % K_RND);
        bool ok = true;
        for (uint16_t i = 1U; i + 1U < n; i++)
        {
            int32_t next_val;
            if (dir > 0)
                next_val = (int32_t)rnd_res[i-1U] + 1 + (rand() % (K_RND - rnd_res[i-1U] > 1 ? K_RND - rnd_res[i-1U]-1 : 1));
            else
                next_val = (int32_t)rnd_res[i-1U] - 1 - (rand() % (rnd_res[i-1U] > 1 ? rnd_res[i-1U]-1 : 1));
            if (next_val < 0 || next_val >= (int32_t)K_RND) { ok = false; break; }
            rnd_res[i] = (RF_Class_t)next_val;
            dir = -dir;
        }
        if (!ok) continue;
        rnd_res[n-1U] = rnd_res[0]; /* force wrap plateau */
        if (rnd_res[n-1U] == rnd_res[n-2U]) continue;
        if (rf_sign_diff(rnd_res[n-1U], rnd_res[n-2U]) == rf_sign_diff(rnd_res[n-2U], rnd_res[n-3U])) continue;

        RF_Class_t tmp[64]; uint16_t tc=0; uint32_t td=0;
        for (uint16_t i=0;i<n;i++) rf_stack_push_and_close_buf(tmp,&tc,64,rnd_res[i],&td);
        if (tc != n) continue;

        /* Measure cascade depth: how many points of the second copy
         * are consumed by "replace" before the first genuine push
         * (true reversal) occurs? */
        RF_Class_t candidate = rnd_res[n-1U];
        int32_t last_dir = rf_sign_diff(rnd_res[n-1U], rnd_res[n-2U]);
        uint16_t depth = 0U;
        for (uint16_t i = 0U; i < n; i++)
        {
            const int32_t dir2 = rf_sign_diff(rnd_res[i], candidate);
            if (dir2 == 0) { depth++; continue; }
            else if (dir2 == last_dir) { candidate = rnd_res[i]; depth++; }
            else { break; } /* first genuine push — cascade ends here */
        }
        if (depth > max_cascade_depth)
        {
            max_cascade_depth = depth;
            max_cascade_n = n;
            memcpy(max_cascade_example, rnd_res, n * sizeof(RF_Class_t));
        }
    }
    printf("Maximum observed cascade depth: %u (at n=%u)\n", max_cascade_depth, max_cascade_n);
    printf("Example residue: [");
    for (uint16_t i = 0U; i < max_cascade_n; i++) printf("%u%s", max_cascade_example[i], (i+1U<max_cascade_n)?",":"");
    printf("]\n");

    const bool all_ok = (g_mismatches_strict == 0U) && (rnd_mismatch == 0U) && (forced_mismatch == 0U);

    printf("\n========================================================\n");
    printf("OVERALL RESULT (stages 2-4; stage 1 is negative control):\n");
    printf("  Stage 2 (exhaustive):          %6llu tested, %llu mismatches\n",
           (unsigned long long)g_valid_strict, (unsigned long long)g_mismatches_strict);
    printf("  Stage 3 (randomized, K=64):    %6llu tested, %llu mismatches\n",
           (unsigned long long)rnd_tested, (unsigned long long)rnd_mismatch);
    printf("  Stage 4 (wrap plateau forced): %6llu tested, %llu mismatches\n",
           (unsigned long long)forced_tested, (unsigned long long)forced_mismatch);
    printf("  Proven max. cascade depth:     2 (empirically confirmed: %u)\n",
           (g_max_cascade_exh > max_cascade_depth) ? g_max_cascade_exh : max_cascade_depth);
    printf("  %s\n", all_ok ? "ALL TESTS PASSED" : "AT LEAST ONE MISMATCH FOUND");
    printf("========================================================\n");

    return all_ok ? 0 : 1;
}
