/*
 *   |     .-.
 *   |    /   \         .-.
 *   |   /     \       /   \       .-.     .-.     _   _
 *   +--/-------\-----/-----\-----/---\---/---\---/-\-/-\/\/---
 *   | /         \   /       \   /     '-'     '-'
 *   |/           '-'         '-'
 *
 *          ____  ___    _____   __________    ____ _       __
 *         / __ \/   |  /  _/ | / / ____/ /   / __ \ |     / /
 *        / /_/ / /| |  / //  |/ / /_  / /   / / / / | /| / / 
 *       / _, _/ ___ |_/ // /|  / __/ / /___/ /_/ /| |/ |/ /  
 *      /_/ |_/_/  |_/___/_/ |_/_/   /_____/\____/ |__/|__/   
 *
 *    Rainflow Counting Algorithm (4-point-method), C99 compliant
 * 
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
 * This implementation uses the 4-point algorithm mentioned in [3,4] and the 3-point HCM method proposed in [2].
 * To take the residue into account, you may implement a custom method or use some
 * predefined functions.
 * 
 * References:
 * [1] ASTM Standard E 1049, 1985 (2011). 
 *     "Standard Practices for Cycle Counting in Fatigue Analysis."
 *     West Conshohocken, PA: ASTM International, 2011.
 * [2] Rainflow - HCM
 *     "Ein Hysteresisschleifen-Zaehlalgorithmus auf werkstoffmechanischer Grundlage"
 *     U.H. Clormann, T. Seeger
 *     1985 TU Darmstadt, Fachgebiet Werkstoffmechanik
 * [3] FVA-Richtlinie, 2010.
 *     "Zaehlverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"
 *     [https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf]
 * [4] Siemens Product Lifecycle Management Software Inc., 2018. 
 *     [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
 * [5] G.Marsh on: "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"
 *     International Journal of Fatigue 82 (2016) 757-765,
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
#include <limits.h>  /* ULLONG_MAX */
#include <stddef.h>  /* size_t, NULL */
#include <stdlib.h>  /* calloc(), free(), abs() */
#include "config.h"  /* Configuration */


#ifndef RFC_VALUE_TYPE
#define RFC_VALUE_TYPE double
#endif

#ifdef RFC_USE_INTEGRAL_COUNTS
#define RFC_COUNTS_VALUE_TYPE    unsigned long long
#define RFC_FULL_CYCLE_INCREMENT (2)
#define RFC_HALF_CYCLE_INCREMENT (1)
#define RFC_COUNTS_LIMIT         (ULLONG_MAX - RFC_FULL_CYCLE_INCREMENT) /* ~18e18 (eff. ~9e18)*/
#else
#define RFC_COUNTS_VALUE_TYPE    double
#define RFC_FULL_CYCLE_INCREMENT (1.0)
#define RFC_HALF_CYCLE_INCREMENT (0.5)
#define RFC_COUNTS_LIMIT         (4.5e15 - RFC_FULL_CYCLE_INCREMENT)
#endif

#ifndef RFC_USE_DELEGATES
#define RFC_USE_DELEGATES 0
#endif

#ifndef RFC_HCM_SUPPORT
#define RFC_HCM_SUPPORT 0
#endif

#ifndef RFC_TP_SUPPORT
#define RFC_TP_SUPPORT 1
#endif

#ifndef RFC_DH_SUPPORT
#define RFC_DH_SUPPORT 0
#endif

#ifndef RFC_GLOBAL_EXTREMA
#define RFC_GLOBAL_EXTREMA 0
#endif


/* Memory allocation functions typedef */
typedef void * ( *rfc_mem_alloc_fcn_t )( void *, size_t num, size_t size );

/* Typedefs */
typedef RFC_VALUE_TYPE          RFC_value_type;      /** Input data value type */
typedef RFC_COUNTS_VALUE_TYPE   RFC_counts_type;     /** Type of counting values */
typedef struct rfc_ctx          rfc_ctx_s;           /** Forward declaration (rainflow context) */
typedef struct rfc_value_tuple  rfc_value_tuple_s;   /** Tuple of value and index position */
typedef struct rfc_class_param  rfc_class_param_s;   /** Class parameters (width, offset, count) */


/* Core functions */
#if !RFC_TP_SUPPORT
bool RFC_init                 ( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                                           RFC_value_type hysteresis );
#else /*RFC_TP_SUPPORT*/
bool RFC_init                 ( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                                           RFC_value_type hysteresis,
                                           rfc_value_tuple_s *tp, size_t tp_cap );
#endif /*RFC_TP_SUPPORT*/
void RFC_deinit               ( void *ctx );
bool RFC_feed                 ( void *ctx, const RFC_value_type* data, size_t count );
bool RFC_feed_tuple           ( void *ctx, rfc_value_tuple_s *data, size_t count );
bool RFC_finalize             ( void *ctx, int residual_method );

#if RFC_USE_DELEGATES
/* Delegates typedef */
typedef  void                ( *rfc_cycle_find_fcn_t )    ( rfc_ctx_s * );
typedef  double              ( *rfc_damage_calc_fcn_t )   ( rfc_ctx_s *, unsigned from_class, unsigned to_class );
typedef  bool                ( *rfc_finalize_fcn_t )      ( rfc_ctx_s *, int residual_methods );
#if RFC_TP_SUPPORT
typedef  rfc_value_tuple_s * ( *rfc_tp_next_fcn_t )       ( rfc_ctx_s *, const rfc_value_tuple_s * );
typedef  bool                ( *rfc_tp_add_fcn_t )        ( rfc_ctx_s *, rfc_value_tuple_s * );
#endif /*RFC_TP_SUPPORT*/
#endif /*RFC_USE_DELEGATES*/


/* Value info struct */
typedef struct rfc_value_tuple
{
    RFC_value_type                      value;                      /**< Value. Don't change order, value field must be first! */
    unsigned                            class;                      /**< Class number, base 0 */
    size_t                              pos;                        /**< Absolute position in input data stream, base 1 */
#if RFC_TP_SUPPORT
    size_t                              tp_pos;                     /**< Position in tp storage, base 0 */
    double                              damage;                     /**< Damage accumulated to this turning point */
#endif /*RFC_TP_SUPPORT*/
} rfc_value_tuple_s;

typedef struct rfc_class_param
{
    unsigned                            count;                      /**< Class count */
    RFC_value_type                      width;                      /**< Class width */
    RFC_value_type                      offset;                     /**< Lower bound of first class */
} rfc_class_param_s;


/**
 * Rainflow context (ctx)
 */
typedef struct rfc_ctx
{
    size_t                              version;                    /**< Version number as sizeof(struct rfctx..), must be 1st field! */

    enum
    {
        RFC_STATE_INIT0,                                            /**< Initialized with zeros */
        RFC_STATE_INIT,                                             /**< Initialized, memory allocated */
        RFC_STATE_BUSY,                                             /**< In counting state */
        RFC_STATE_BUSY_INTERIM,                                     /**< In counting state, having still one interim turning point (not included) */
        RFC_STATE_FINALIZE,                                         /**< Finalizing */
        RFC_STATE_FINISHED,                                         /**< Counting finished, memory still allocated */
        RFC_STATE_ERROR,                                            /**< An error occured */
    }                                   state;                      /**< Current counting state */

    enum
    {
        RFC_ERROR_INVARG,
        RFC_ERROR_MEMORY,
    }                                   error;                      /**< Error code */

    enum
    {
        RFC_FLAGS_COUNT_MATRIX          = 1 << 0,                   /**< Count into matrix */
        RFC_FLAGS_COUNT_RP              = 1 << 1,                   /**< Count into range pair */
        RFC_FLAGS_COUNT_LC_UP           = 1 << 2,                   /**< Count into level crossing (only rising slopes) */
        RFC_FLAGS_COUNT_LC_DN           = 1 << 3,                   /**< Count into level crossing (only falling slopes) */
        RFC_FLAGS_COUNT_LC              = RFC_FLAGS_COUNT_LC_UP     /**< Count into level crossing (all slopes) */
                                        | RFC_FLAGS_COUNT_LC_DN,
        RFC_FLAGS_COUNT_ALL             = RFC_FLAGS_COUNT_MATRIX    /**< Count all */
                                        | RFC_FLAGS_COUNT_RP
                                        | RFC_FLAGS_COUNT_LC,
        RFC_FLAGS_ENFORCE_MARGIN        = 1 << 8,                   /**< Enforce first and last data point are turning points */
    }
                                        flags;                      /**< Flags */
    enum
    {
        RFC_COUNTING_METHOD_UNKNOWN     = -1,                       /**< Method is unknown */
        RFC_COUNTING_METHOD_NONE        =  0,                       /**< No counting */
        RFC_COUNTING_METHOD_4PTM        =  1,                       /**< 4 point algorithm (default) */
#if RFC_HCM_SUPPORT
        RFC_COUNTING_METHOD_HCM         =  2,                       /**< 3 point algorithm, Clormann/Seeger (HCM) method */
#endif /*RFC_HCM_SUPPORT*/
    }
                                        counting_method;
    enum 
    {
        RFC_RES_NONE                    = 0,                        /**< No residual method */
        RFC_RES_IGNORE,                                             /**< Ignore residue (same as RFC_RES_NONE) */
        RFC_RES_DISCARD,                                            /**< Discard residue (empty residue) */
        RFC_RES_HALFCYCLES,                                         /**< ASTM */
        RFC_RES_FULLCYCLES,                                         /**< Count half cycles as full cycles */
        RFC_RES_CLORMANN_SEEGER,                                    /**< Clormann/Seeger method */
        RFC_RES_REPEATED,                                           /**< Repeat residue and count closed cycles */
        RFC_RES_RP_DIN45667,                                        /**< Count residue according to range pair in DIN-45667 */
    }
                                        residual_method;
#if RFC_DH_SUPPORT
    enum
    {
        RFC_SD_NONE                     = -1,                       /**< No spread damage calculation */
        RFC_SD_HALF_23                  =  0,                       /**< Equally split damage between P2 and P3 */
        RFC_SD_RAMP_AMPLITUDE_23        =  1,                       /**< Spread damage according to amplitude over points between P2 and P3 */
        RFC_SD_RAMP_DAMAGE_23           =  2,                       /**< Spread damage linearly over points between P2 and P3 */
        RFC_SD_RAMP_AMPLITUDE_24        =  3,                       /**< Spread damage according to amplitude over points between P2 and P4 */  
        RFC_SD_RAMP_DAMAGE_24           =  4,                       /**< Spread damage linearly over points between P2 and P4 */
        RFC_SD_FULL_P2                  =  5,                       /**< Assign damage to P2 */
        RFC_SD_FULL_P3                  =  6,                       /**< Assign damage to P3 */
        RFC_SD_TRANSIENT_23             =  7,                       /**< Spread damage transient according to amplitude over points between P2 and P3 */
        RFC_SD_TRANSIENT_23c            =  8,                       /**< Spread damage transient according to amplitude over points between P2 and P4 only until cycle is closed */
    }
                                        spread_damage_method;
#endif /*RFC_DH_SUPPORT*/

    /* Memory allocation functions */
    rfc_mem_alloc_fcn_t                 mem_alloc;                  /**< Allocate initialized memory */

    /* Counter increments */
    RFC_counts_type                     full_inc;                   /**< Increment for a full cycle */
    RFC_counts_type                     half_inc;                   /**< Increment for a half cycle, used by some residual algorithms */
    RFC_counts_type                     curr_inc;                   /**< Current increment, used by counting algorithms */

    /* Rainflow class parameters */
    unsigned                            class_count;                /**< Class count */
    RFC_value_type                      class_width;                /**< Class width */
    RFC_value_type                      class_offset;               /**< Lower bound of first class */
    RFC_value_type                      hysteresis;                 /**< Hysteresis filtering */

    /* Woehler curve */
    double                              wl_sd;                      /**< Fatigue resistance range (amplitude) */
    double                              wl_nd;                      /**< Cycles at wl_sd */
    double                              wl_k;                       /**< Woehler gradient above wl_sd */
    double                              wl_k2;                      /**< Woehler gradient below wl_sd */
    double                              wl_omission;                /**< Omission level */

#if RFC_USE_DELEGATES
    /* Delegates (optional, may be NULL) */
#if RFC_TP_SUPPORT
    rfc_tp_next_fcn_t                   tp_next_fcn;                /**< Test if new turning point exists */
    rfc_tp_add_fcn_t                    tp_add_fcn;                 /**< Handling new turning points */
#endif /*RFC_TP_SUPPORT*/
    rfc_finalize_fcn_t                  finalize_fcn;               /**< Finalizing function */
    rfc_cycle_find_fcn_t                cycle_find_fcn;             /**< Find next cycle(s) and process */
    rfc_damage_calc_fcn_t               damage_calc_fcn;            /**< Damage calculating function */
#endif /*RFC_USE_DELEGATES*/
    
    /* Residue */
    rfc_value_tuple_s                  *residue;                    /**< Buffer for residue */
    size_t                              residue_cap;                /**< Buffer capacity in number of elements (max. 2*class_count) */
    size_t                              residue_cnt;                /**< Number of value tuples in buffer */

    /* Non-sparse storages (optional, may be NULL) */
    RFC_counts_type                    *matrix;                     /**< Rainflow matrix */
    RFC_counts_type                    *rp;                         /**< Range pair counts */
    RFC_counts_type                    *lc;                         /**< Level crossing counts */

#if RFC_TP_SUPPORT
    /* Turning points storage (optional, may be NULL) */
    rfc_value_tuple_s                  *tp;                         /**< Buffer for turning points, pointer may be changed thru memory reallocation! */
    size_t                              tp_cap;                     /**< Buffer capacity (number of elements) */
    size_t                              tp_cnt;                     /**< Number of turning points in buffer */
    bool                                tp_locked;                  /**< If tp_locked, tp is freezed */
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    double                             *dh;                         /**< Damage history, pointer may be changed thru memory reallocation! */
    size_t                              dh_cap;                     /**< Capacity of dh */
    size_t                              dh_cnt;                     /**< Number of values in dh */
#endif /*RFC_DH_SUPPORT*/

    /* Damage */
    double                             *damage_lut;                 /**< Damage look-up table */
    double                              pseudo_damage;              /**< Cumulated pseudo damage */
    
    /* Internal usage */
    struct internal
    {
        int                             slope;                      /**< Current signal slope */
        rfc_value_tuple_s               extrema[2];                 /**< Local or global extrema depending on RFC_GLOBAL_EXTREMA */
        bool                            extrema_changed;            /**< True if one extrema has changed */
        size_t                          pos;                        /**< Absolute position in data input stream, base 1 */
#if RFC_TP_SUPPORT
        rfc_value_tuple_s               margin[2];                  /**< First and last data point */
        rfc_value_tuple_s               tp_delayed;                 /**< Delay stage when RFC_FLAGS_ENFORCE_MARGIN is set */
#endif /*RFC_TP_SUPPORT*/
#if RFC_HCM_SUPPORT
        struct hcm
        {
            /* Residue */
            rfc_value_tuple_s          *stack;                      /**< Stack */
            size_t                      stack_cap;                  /**< Stack capacity in number of elements (max. 2*class_count) */
            int                         IR;                         /**< Pointer to residue stack, first turning point of cycles able to close, base 1 */
            int                         IZ;                         /**< Pointer to residue stack, last turning point of cycles able to close, base 1 */
        }                               hcm;
#endif /*RFC_HCM_SUPPORT*/
    }
                                        internal;
} rfc_ctx_s;
