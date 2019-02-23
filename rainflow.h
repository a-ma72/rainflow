/*
 *
 *   |                     .-.
 *   |                    /   \
 *   |     .-.===========/     \         .-.
 *   |    /   \         /       \       /   \
 *   |   /     \       /         \     /     \         .-.
 *   +--/-------\-----/-----------\---/-------\-------/---\
 *   | /         \   /             '-'=========\     /     \   /
 *   |/           '-'                           \   /       '-'
 *   |                                           '-'
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
 * [1] "Standard Practices for Cycle Counting in Fatigue Analysis."
 *     ASTM Standard E 1049, 1985 (2011). 
 *     West Conshohocken, PA: ASTM International, 2011.
 * [2] "Rainflow - HCM / Ein Hysteresisschleifen-Zaehlalgorithmus auf werkstoffmechanischer Grundlage"
 *     U.H. Clormann, T. Seeger
 *     1985 TU Darmstadt, Fachgebiet Werkstoffmechanik
 * [3] "Zaehlverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"
 *     FVA-Richtlinie, 2010.
 *     [https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf]
 * [4] Siemens Product Lifecycle Management Software Inc., 2018. 
 *     [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
 * [5] "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"
 *     G.Marsh;
 *     International Journal of Fatigue 82 (2016) 757-765,
 *     [https://doi.org/10.1016/j.ijfatigue.2015.10.007]
 * [6] "Betriebsfestigkeit - Verfahren und Daten zur Bauteilberechnung"
 *     Haibach, Erwin; Springer Verlag
 * []  "Schaedigungsbasierte Hysteresefilter"; Hack, M, D386 (Diss Univ. Kaiserslautern), Shaker Verlag Aachen, 1998, ISBN 3-8265-3936-2
 * []  "Hysteresis and Phase Transition"
 *     Brokate, M.; Sprekels, J.; Applied Mathematical Sciences 121; Springer, New York, 1996
 * []  "Rainflow counting and energy dissipation in elastoplasticity"; Eur. J. Mech. A/Solids 15, pp. 705-737, 1996
 *     Brokate, M.; Dressler, K.; Krejci, P.
 * []  "Multivariate Density Estimation: Theory, Practice and Visualization". New York, Chichester, Wiley & Sons, 1992
 *     Scott, D.
 * []  "Werkstoffmechanik - Bauteile sicher beurteilen undWerkstoffe richtig einsetzen"; 
 *      Ralf Buergel, Hans Albert Richard, Andre Riemer; Springer FachmedienWiesbaden 2005, 2014
 * []  "Zaelverfahren und Lastannahme in der Betriebsfestigkeit";
 *     Michael Koehler, Sven Jenne / Kurt Poetter, Harald Zenner; Springer-Verlag Berlin Heidelberg 2012
 *
 *
 *================================================================================
 * BSD 2-Clause License
 * 
 * Copyright (c) 2019, Andras Martin
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

#ifndef RAINFLOW_H
#define RAINFLOW_H

/* This version is generated via coan (http://coan2.sourceforge.net/) */

#include "config.h"  /* Configuration */

#ifndef RFC_VALUE_TYPE
#define RFC_VALUE_TYPE double
#endif /*RFC_VALUE_TYPE*/

#ifdef RFC_USE_INTEGRAL_COUNTS
#define RFC_COUNTS_VALUE_TYPE    unsigned long long
#define RFC_FULL_CYCLE_INCREMENT (2)
#define RFC_HALF_CYCLE_INCREMENT (1)
#define RFC_COUNTS_LIMIT         (ULLONG_MAX - RFC_FULL_CYCLE_INCREMENT) /* ~18e18 (eff. ~9e18)*/
#else /*RFC_USE_INTEGRAL_COUNTS*/
#define RFC_COUNTS_VALUE_TYPE    double
#define RFC_FULL_CYCLE_INCREMENT (1.0)
#define RFC_HALF_CYCLE_INCREMENT (0.5)
#define RFC_COUNTS_LIMIT         (4.5e15 - RFC_FULL_CYCLE_INCREMENT)
#endif /*RFC_USE_INTEGRAL_COUNTS*/



// Notes on mix C and C++ headers:
// https://developers.redhat.com/blog/2016/02/29/why-cstdlib-is-more-complicated-than-you-might-think/
// Avoid including C standard headers in a C++ namespace!
#ifdef __cplusplus
#include <cstdbool>  /* bool, true, false */
#include <cstdint>   /* ULLONG_MAX */
#include <climits>   /* ULLONG_MAX */
#include <cstddef>   /* size_t, NULL */
#ifndef RFC_CPP_NAMESPACE
#define RFC_CPP_NAMESPACE rainflow_C
#endif /*RFC_CPP_NAMESPACE*/
namespace RFC_CPP_NAMESPACE {
#else /*!__cplusplus*/
#include <stdbool.h> /* bool, true, false */
#include <stdint.h>  /* ULLONG_MAX */
#include <limits.h>  /* ULLONG_MAX */
#include <stddef.h>  /* size_t, NULL */
#endif /*__cplusplus*/


/* Memory allocation aim info */
enum rfc_mem_aim
{
    RFC_MEM_AIM_TEMP                =  0,                           /**< Error on accessing memory for temporary storage */
    RFC_MEM_AIM_RESIDUE             =  1,                           /**< Error on accessing memory for residue */
    RFC_MEM_AIM_MATRIX              =  2,                           /**< Error on accessing memory for rf matrix */
    RFC_MEM_AIM_RP                  =  3,                           /**< Error on accessing memory for range pair counting */
    RFC_MEM_AIM_LC                  =  4,                           /**< Error on accessing memory for level crossing */
};


/* Flags */
enum rfc_flags
{
    RFC_FLAGS_DEFAULT               = -1,
    RFC_FLAGS_COUNT_RFM             =  1 << 0,                      /**< Count into rainflow matrix */
    RFC_FLAGS_COUNT_DAMAGE          =  1 << 1,                      /**< Count damage */
    RFC_FLAGS_COUNT_ALL             =  RFC_FLAGS_COUNT_RFM          /**< Count all */
                                    |  RFC_FLAGS_COUNT_DAMAGE
};


enum rfc_debug_flags
{
    RFC_FLAGS_LOG_CLOSED_CYCLES     =  1 << 0,                      /**< Log closed cycles */
};


enum rfc_state
{
    RFC_STATE_INIT0,                                                /**< Initialized with zeros */
    RFC_STATE_INIT,                                                 /**< Initialized, memory allocated */
    RFC_STATE_BUSY,                                                 /**< In counting state */
    RFC_STATE_BUSY_INTERIM,                                         /**< In counting state, having still one interim turning point (not included) */
    RFC_STATE_FINALIZE,                                             /**< Finalizing */
    RFC_STATE_FINISHED,                                             /**< Counting finished, memory still allocated */
    RFC_STATE_ERROR,                                                /**< An error occurred */
};


enum rfc_error
{
    RFC_ERROR_UNEXP                 = -1,                           /**< Unexpected error */
    RFC_ERROR_NOERROR               =  0,                           /**< No error */
    RFC_ERROR_INVARG                =  1,                           /**< Invalid arguments passed */
    RFC_ERROR_UNSUPPORTED           =  2,                           /**< Unsupported feature */
    RFC_ERROR_MEMORY                =  3,                           /**< Error on memory allocation */
};


enum rfc_res_method
{
    /* Don't change order! */
    RFC_RES_NONE                    = 0,                            /**< No residual method */
    RFC_RES_IGNORE                  = 1,                            /**< Ignore residue (same as RFC_RES_NONE) */
    RFC_RES_COUNT                                                   /**< Number of options */
};


/* Typedefs */
typedef                 RFC_VALUE_TYPE          rfc_value_t;                /** Input data value type */
typedef                 RFC_COUNTS_VALUE_TYPE   rfc_counts_t;               /** Type of counting values */
typedef     struct      rfc_value_tuple         rfc_value_tuple_s;          /** Tuple of value and index position */
typedef     struct      rfc_ctx                 rfc_ctx_s;                  /** Forward declaration (rainflow context) */
typedef     enum        rfc_mem_aim             rfc_mem_aim_e;              /** Memory accessing mode */
typedef     enum        rfc_flags               rfc_flags_e;                /** Flags, see RFC_FLAGS... */
typedef     enum        rfc_debug_flags         rfc_debug_flags_e;          /** Flags, see RFC_DEBUG_FLAGS... */
typedef     enum        rfc_state               rfc_state_e;                /** Counting state, see RFC_STATE... */
typedef     enum        rfc_error               rfc_error_e;                /** Recent error, see RFC_ERROR... */
typedef     enum        rfc_res_method          rfc_res_method_e;           /** Method when count residue into matrix, see RFC_RES... */

/* Memory allocation functions typedef */
typedef     void *   ( *rfc_mem_alloc_fcn_t )   ( void *, size_t num, size_t size, rfc_mem_aim_e aim );     /** Memory allocation functor */

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

/* Core functions */
bool    RFC_init                    (       void *ctx, unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                                                       rfc_value_t hysteresis, rfc_flags_e flags );
bool    RFC_wl_init_elementary      (       void *ctx, double sx, double nx, double k );
bool    RFC_deinit                  (       void *ctx );
bool    RFC_feed                    (       void *ctx, const rfc_value_t* data, size_t count );
bool    RFC_finalize                (       void *ctx, rfc_res_method_e residual_method );


#ifdef __cplusplus
}  // extern "C"
#endif /*__cplusplus*/


/* Value info struct */
struct rfc_value_tuple
{
    rfc_value_t                         value;                      /**< Value. Don't change order, value field must be first! */
    unsigned                            cls;                        /**< Class number, base 0 */
    size_t                              pos;                        /**< Absolute position in input data stream, base 1 */
};


/**
 * Rainflow context (ctx)
 */
struct rfc_ctx
{
    size_t                              version;                    /**< Version number as sizeof(struct rfctx..), must be 1st field! */

    /* State and error information */
    rfc_state_e                         state;                      /**< Current counting state */
    rfc_error_e                         error;                      /**< Error code */

    /* Methods */
    rfc_res_method_e                    residual_method;            /**< Used on finalizing */

    /* Memory allocation functions */
    rfc_mem_alloc_fcn_t                 mem_alloc;                  /**< Allocate initialized memory */

    /* Counter increments */
    rfc_counts_t                        full_inc;                   /**< Increment for a full cycle */
    rfc_counts_t                        half_inc;                   /**< Increment for a half cycle, used by some residual algorithms */
    rfc_counts_t                        curr_inc;                   /**< Current increment, used by counting algorithms */

    /* Rainflow class parameters */
    unsigned                            class_count;                /**< Class count */
    rfc_value_t                         class_width;                /**< Class width */
    rfc_value_t                         class_offset;               /**< Lower bound of first class */
    rfc_value_t                         hysteresis;                 /**< Hysteresis filtering, slope must exceed hysteresis to be counted! */

    /* Woehler curve */
    double                              wl_sx;                      /**< Sa of any point on the Woehler curve */
    double                              wl_nx;                      /**< Cycles for Sa on the Woehler curve */
    double                              wl_k;                       /**< Woehler slope, always negative */

    /* Residue */
    rfc_value_tuple_s                  *residue;                    /**< Buffer for residue */
    size_t                              residue_cap;                /**< Buffer capacity in number of elements (max. 2*class_count) */
    size_t                              residue_cnt;                /**< Number of value tuples in buffer */

    /* Non-sparse storages (optional, may be NULL) */
    rfc_counts_t                       *rfm;                        /**< Rainflow matrix, always class_count^2 elements (row-major, row=from, to=col). */

    /* Damage */
    double                              damage;                     /**< Cumulated damage */
    double                              damage_residue;             /**< Partial damage in .damage influenced by taking residue into account (after finalizing) */

    /* Internal usage */
    struct internal
    {
        int                             flags;                      /**< Flags (enum rfc_flags) */
        int                             slope;                      /**< Current signal slope */
        rfc_value_tuple_s               extrema[2];                 /**< Local or global extrema depending on RFC_GLOBAL_EXTREMA */
        size_t                          pos;                        /**< Absolute position in data input stream, base 1 */
        size_t                          pos_offset;                 /**< Offset for pos */
        rfc_value_tuple_s               residue[3];                 /**< Static residue (if class_count is zero) */
        size_t                          residue_cap;                /**< Capacity of static residue */
        bool                            res_static;                 /**< true, if .residue refers the static residue .internal.residue */
    }
                                        internal;
};

#ifdef __cplusplus
}  // namespace RFC_CPP_NAMESPACE
#endif /*__cplusplus*/

#endif /*RAINFLOW_H*/