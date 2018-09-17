/** Rainflow Counting Algorithm (4-point-method), C89 compliant
 * 
 * "Rainflow Counting" consists of four main steps:
 *    1. Hysteresis Filtering
 *    2. Peak-Valley Filtering
 *    3. Discretization
 *    4. Four Point Counting Method:
 *    
 *                     * D
 *                    / \       Closed, if min(B,C) >= min(A,D) && max(B,C) <= max(A,D)
 *             B *<--/          Slope B-C is counted and removed from residue
 *              / \ /
 *             /   * C
 *          \ /
 *           * A
 *
 * These steps are fully documented in standards such as 
 * ASTM E1049 "Standard Practices for Cycle Counting in Fatigue Analysis" [1].
 * This implementation uses the 4-point algorithm mentioned in [2].
 * To take the residue into account, you may implement a custom method or use some
 * predefined functions.
 * 
 * References:
 * [1] ASTM Standard E 1049, 1985 (2011). 
 *     "Standard Practices for Cycle Counting in Fatigue Analysis."
 *     West Conshohocken, PA: ASTM International, 2011.
 * [2] [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
 * [3] G.Marsh on: "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"
 *     International Journal of Fatigue 82 (2016) 757–765,
 *     [https://doi.org/10.1016/j.ijfatigue.2015.10.007]
 * []  Hack, M: Schaedigungsbasierte Hysteresefilter; D386 (Diss Univ. Kaiserslautern), Shaker Verlag Aachen, 1998, ISBN 3-8265-3936-2
 * []  Brokate, M; Sprekels, J, Hysteresis and Phase Transition, Applied Mathematical Sciences 121, Springer,  New York, 1996
 * []  Brokate, M; Dressler, K; Krejci, P: Rainflow counting and energy dissipation in elastoplasticity, Eur. J. Mech. A/Solids 15, pp. 705-737, 1996
 * []  Scott, D.: Multivariate Density Estimation: Theory, Practice and Visualization. New York, Chichester, Wiley & Sons, 1992
 *
 *                                                                                                                                                          *
 *================================================================================
 * BSD 2-Clause License
 * 
 * Copyright (c) 2018, Andras Martin
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *================================================================================
 */


#include <stdbool.h> /* bool, true, false */
#include <stdint.h>  /* ULLONG_MAX */
#include <stddef.h>  /* size_t, NULL */
#include <stdlib.h>  /* calloc(), free(), abs() */


#ifndef RFC_VALUE_TYPE
#define RFC_VALUE_TYPE double
#endif

#ifdef RFC_USE_INTEGRAL_COUNTS
#define RFC_COUNTS_VALUE_TYPE    (unsigned long long)
#define RFC_FULL_CYCLE_INCREMENT (2)
#define RFC_HALF_CYCLE_INCREMENT (1)
#define RFC_COUNTS_LIMIT         (ULLONG_MAX - RFC_FULL_CYCLE_INCREMENT) /* ~18e18 (eff. ~9e18)*/
#else
#define RFC_COUNTS_VALUE_TYPE    double
#define RFC_FULL_CYCLE_INCREMENT (1.0)
#define RFC_HALF_CYCLE_INCREMENT (0.5)
#define RFC_COUNTS_LIMIT         (4.5e15 - RFC_FULL_CYCLE_INCREMENT)
#endif

typedef RFC_VALUE_TYPE          RFC_value_type;  /** Input data value type */
typedef RFC_COUNTS_VALUE_TYPE   RFC_counts_type; /** Type of counting values */
typedef struct rfctx            rfctx_s;         /** Forward declaration (rainflow context) */
typedef struct value_tuple      value_tuple_s;   /** Tuple of value and index position */


/* Core */
bool RFC_init                 ( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                                           RFC_value_type hysteresis );
void RFC_feed                 ( void *ctx, const RFC_value_type* data, size_t count, bool do_finalize );
void RFC_finalize             ( void *ctx );
void RFC_deinit               ( void *ctx );

/* Memory allocation functions typedef */
typedef void * ( *rfc_mem_calloc_fcn_t )    ( size_t num, size_t size );
typedef void   ( *rfc_mem_free_fcn_t )      ( void * );

/* Value info struct */
typedef struct value_tuple
{
    RFC_value_type                  value;                      /**< Value. Don't change order, value field must be first! */
    unsigned                        class;                      /**< Class number, base 0 */
} value_tuple_s;


/**
 * Rainflow context
 */
typedef struct rfctx
{
    size_t                          version;                    /**< Version number as sizeof(struct rfctx..), must be 1st field! */

    enum
    {
        RFC_STATE_INIT0,                                        /**< Initialized with zeros */
        RFC_STATE_INIT,                                         /**< Initialized, memory allocated */
        RFC_STATE_BUSY,                                         /**< In counting state */
        RFC_STATE_FINALIZE,                                     /**< Finalizing */
        RFC_STATE_FINISHED,                                     /**< Counting finished, memory still allocated */
        RFC_STATE_ERROR,                                        /**< An error occured */
    }                               state;                      /**< Current counting state */

    enum
    {
        RFC_FLAGS_COUNT_MATRIX      = 1 << 0,                   /**< Count into matrix */
        RFC_FLAGS_COUNT_RP          = 1 << 1,                   /**< Count into range pair */
        RFC_FLAGS_COUNT_LC_UP       = 1 << 2,                   /**< Count into level crossing (only rising slopes) */
        RFC_FLAGS_COUNT_LC_DN       = 1 << 3,                   /**< Count into level crossing (only falling slopes) */
        RFC_FLAGS_COUNT_LC          = RFC_FLAGS_COUNT_LC_UP     /**< Count into level crossing (all slopes) */
                                    | RFC_FLAGS_COUNT_LC_DN,
        RFC_FLAGS_COUNT_ALL         = -1                        /** Count all */
    }
                                    flags;                      /**< Flags */

    /* Counter increments */
    RFC_counts_type                 full_inc;                   /**< Increment for a full cycle */
    RFC_counts_type                 half_inc;                   /**< Increment for a half cycle, used by some residual algorithms */
    RFC_counts_type                 curr_inc;                   /**< Current increment, used by counting algorithms */

    /* Rainflow class parameters */
    unsigned                        class_count;                /**< Class count */
    RFC_value_type                  class_width;                /**< Class width */
    RFC_value_type                  class_offset;               /**< Lower bound of first class */
    RFC_value_type                  hysteresis;                 /**< Hysteresis Filtering */

    /* Woehler curve */
    double                          wl_sd;                      /**< Fatigue resistance range */
    double                          wl_nd;                      /**< Cycles at wl_sd */
    double                          wl_k;                       /**< Woehler gradient */

    /* Memory allocation functions */
    rfc_mem_calloc_fcn_t            mem_alloc;                  /**< Allocate initialized memory */
    rfc_mem_free_fcn_t              mem_free;                   /**< Free memory, allocated by mem_alloc */

    /* Residue */
    value_tuple_s                  *residue;                    /**< Buffer for residue */
    size_t                          residue_size;               /**< Buffer size (max. 2*class_count) */
    size_t                          residue_cnt;                /**< Buffer content size */

    /* Non-sparse storages (optional, may be NULL) */
    RFC_counts_type                *matrix;                     /**< Rainflow matrix */

    /* Damage */
    double                          pseudo_damage;              /**< Cumulated pseudo damage */
    
    /* Internal usage */
    struct internal
    {
        int                         slope;                      /**< Current signal slope */
        value_tuple_s               extrema[2];                 /**< Local extrema */
        size_t                      pos;                        /**< Absolute position in data input stream, base 1 */
    } internal;
} rfctx_s;
