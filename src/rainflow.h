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
 * Copyright (c) 2022, Andras Martin
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

#if COAN_INVOKED
/* This version is generated via coan (http://coan2.sourceforge.net/) */
#endif /*COAN_INVOKED*/

#if RFC_HAVE_CONFIG_H
#include "config.h"  /* Configuration */
#endif

#ifndef ON
#define ON (1)
#endif /*ON*/
#ifndef OFF
#define OFF (0)
#endif /*OFF*/

#define RFC_CLASS_COUNT_MAX (1024)

#ifndef RFC_VALUE_TYPE
#define RFC_VALUE_TYPE double
#endif /*RFC_VALUE_TYPE*/

#if RFC_USE_INTEGRAL_COUNTS
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

#if RFC_MINIMAL
#undef  RFC_USE_DELEGATES
#define RFC_USE_DELEGATES    OFF
#undef  RFC_HCM_SUPPORT
#define RFC_HCM_SUPPORT      OFF
#undef  RFC_ASTM_SUPPORT
#define RFC_ASTM_SUPPORT     OFF
#undef  RFC_TP_SUPPORT
#define RFC_TP_SUPPORT       OFF
#undef  RFC_DH_SUPPORT
#define RFC_DH_SUPPORT       OFF
#undef  RFC_AT_SUPPORT
#define RFC_AT_SUPPORT       OFF
#undef  RFC_AR_SUPPORT
#define RFC_AR_SUPPORT       OFF
#undef  RFC_GLOBAL_EXTREMA
#define RFC_GLOBAL_EXTREMA   OFF
#undef  RFC_DAMAGE_FAST
#define RFC_DAMAGE_FAST      OFF
#else /*!RFC_MINIMAL*/
#ifndef RFC_MINIMAL
#define RFC_MINIMAL OFF
#endif /*RFC_MINIMAL*/
#ifndef RFC_USE_DELEGATES
#define RFC_USE_DELEGATES OFF
#endif /*RFC_USE_DELEGATES*/
#ifndef RFC_HCM_SUPPORT
#define RFC_HCM_SUPPORT OFF
#endif /*RFC_HCM_SUPPORT*/
#ifndef RFC_ASTM_SUPPORT
#define RFC_ASTM_SUPPORT OFF
#endif /*RFC_ASTM_SUPPORT*/
#ifndef RFC_TP_SUPPORT
#define RFC_TP_SUPPORT ON
#endif /*RFC_TP_SUPPORT*/
#ifndef RFC_DH_SUPPORT
#define RFC_DH_SUPPORT OFF
#endif /*RFC_DH_SUPPORT*/
#ifndef RFC_AT_SUPPORT
#define RFC_AT_SUPPORT OFF
#endif /*RFC_AT_SUPPORT*/
#ifndef RFC_AR_SUPPORT
#define RFC_AR_SUPPORT OFF
#endif /*RFC_AR_SUPPORT*/
#ifndef RFC_GLOBAL_EXTREMA
#define RFC_GLOBAL_EXTREMA OFF
#endif /*RFC_GLOBAL_EXTREMA*/
#ifndef RFC_DAMAGE_FAST
#define RFC_DAMAGE_FAST ON
#endif /*RFC_DAMAGE_FAST*/
#ifndef RFC_DEBUG_FLAGS
#define RFC_DEBUG_FLAGS OFF
#endif /*RFC_DEBUG_FLAGS*/
#endif /*RFC_MINIMAL*/


/* Notes on mix C and C++ headers:
 * https://developers.redhat.com/blog/2016/02/29/why-cstdlib-is-more-complicated-than-you-might-think/
 * Avoid including C standard headers in a C++ namespace! */
#ifdef __cplusplus
#include <cstdbool>  /* bool, true, false */
#include <cstdint>   /* ULLONG_MAX */
#include <climits>   /* ULLONG_MAX */
#include <cstddef>   /* size_t, NULL */
#if RFC_DEBUG_FLAGS
#include <cstdio>
#include <cstdarg>
#endif /*RFC_DEBUG_FLAGS*/
#ifndef RFC_CPP_NAMESPACE
#define RFC_CPP_NAMESPACE rainflow_C
#endif /*RFC_CPP_NAMESPACE*/
namespace RFC_CPP_NAMESPACE {
extern "C" {
#else /*!__cplusplus*/
#include <stdbool.h> /* bool, true, false */
#include <stdint.h>  /* ULLONG_MAX */
#include <limits.h>  /* ULLONG_MAX */
#include <stddef.h>  /* size_t, NULL */
#if RFC_DEBUG_FLAGS
#include <stdio.h>
#include <stdarg.h>
#endif /*RFC_DEBUG_FLAGS*/
#endif /*__cplusplus*/

#pragma pack(push, 1)


/* Memory allocation aim info */
enum rfc_mem_aim
{
    RFC_MEM_AIM_TEMP                =  0,                           /**< Error on accessing memory for temporary storage */
    RFC_MEM_AIM_RESIDUE             =  1,                           /**< Error on accessing memory for residue */
    RFC_MEM_AIM_MATRIX              =  2,                           /**< Error on accessing memory for rf matrix */
    RFC_MEM_AIM_RP                  =  3,                           /**< Error on accessing memory for range pair counting */
    RFC_MEM_AIM_LC                  =  4,                           /**< Error on accessing memory for level crossing */
#if RFC_TP_SUPPORT
    RFC_MEM_AIM_TP                  =  5,                           /**< Error on accessing memory for turning points */
#endif
#if RFC_DAMAGE_FAST
    RFC_MEM_AIM_DLUT                =  6,                           /**< Error on accessing memory for damage look-up table */
#endif
#if RFC_AT_SUPPORT
    RFC_MEM_AIM_ALUT                =  7,                           /**< Error on accessing memory for amplitude look-up table */
#endif
#if RFC_HCM_SUPPORT
    RFC_MEM_AIM_HCM                 =  8,                           /**< Error on accessing memory for HCM algorithm */
#endif
#if RFC_DH_SUPPORT
    RFC_MEM_AIM_DH                  =  9,                           /**< Error on accessing memory for damage history */
#endif
#if !RFC_MINIMAL
    RFC_MEM_AIM_RFM_ELEMENTS        = 10,                           /**< Error on accessing memory for rf matrix elements */
#endif /*!RFC_MINIMAL*/
};


/* Flags */
enum rfc_flags
{
    RFC_FLAGS_DEFAULT               = -1,
    RFC_FLAGS_COUNT_RFM             =  1 << 0,                      /**< Count into rainflow matrix */
    RFC_FLAGS_COUNT_DAMAGE          =  1 << 1,                      /**< Count damage */
#if !RFC_MINIMAL
#if RFC_DH_SUPPORT
    RFC_FLAGS_COUNT_DH              =  1 << 2,                      /**< Spread damage */
#endif /*RFC_DH_SUPPORT*/
    RFC_FLAGS_COUNT_RP              =  1 << 3,                      /**< Count into range pair */
    RFC_FLAGS_COUNT_LC_UP           =  1 << 4,                      /**< Count into level crossing (only rising slopes) */
    RFC_FLAGS_COUNT_LC_DN           =  1 << 5,                      /**< Count into level crossing (only falling slopes) */
    RFC_FLAGS_COUNT_LC              =  RFC_FLAGS_COUNT_LC_UP        /**< Count into level crossing (all slopes) */
                                    |  RFC_FLAGS_COUNT_LC_DN,
    RFC_FLAGS_COUNT_MK              =  1 << 6,                      /**< Live damage counter (Miner consequent) */
    RFC_FLAGS_ENFORCE_MARGIN        =  1 << 7,                      /**< Enforce first and last data point are turning points */
#endif /*!RFC_MINIMAL*/
    RFC_FLAGS_COUNT_ALL             =  RFC_FLAGS_COUNT_RFM          /**< Count all */
                                    |  RFC_FLAGS_COUNT_DAMAGE
#if RFC_DH_SUPPORT
                                    |  RFC_FLAGS_COUNT_DH
#endif /*RFC_DH_SUPPORT*/
#if !RFC_MINIMAL
                                    |  RFC_FLAGS_COUNT_RP
                                    |  RFC_FLAGS_COUNT_LC
                                    |  RFC_FLAGS_COUNT_MK,
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
    RFC_FLAGS_TPPRUNE_PRESERVE_POS  =  1 << 8,                      /**< Preserve stream position information on pruning */
    RFC_FLAGS_TPPRUNE_PRESERVE_RES  =  1 << 9,                      /**< Preserve turning points that exist in resiude on pruning */
    RFC_FLAGS_TPAUTOPRUNE           =  1 << 10,                     /**< Automatic prune on tp */
#endif /*RFC_TP_SUPPORT*/
#if RFC_AR_SUPPORT
    RFC_FLAGS_AUTORESIZE            =  1 << 11,                     /**< Automatically resize buffers for rp, lc, and rfm */
#endif /*RFC_AR_SUPPORT*/
};


enum rfc_debug_flags
{
    RFC_FLAGS_LOG_CLOSED_CYCLES     =  1 << 0,                      /**< Log closed cycles */
#if RFC_TP_SUPPORT
    RFC_FLAGS_LOG_READ_TP           =  1 << 1,                      /**< Read from tp storage */
    RFC_FLAGS_LOG_WRITE_TP          =  1 << 2,                      /**< Write to tp storage */
    RFC_FLAGS_LOG_TP_REFEED         =  1 << 3,                      /**< Trace tp storage while tp_refeed */
#endif /*RFC_TP_SUPPORT*/
};


#if !RFC_MINIMAL
/* See RFC_damage_from_rp() */
enum rfc_rp_damage_method
{
    RFC_RP_DAMAGE_CALC_METHOD_DEFAULT     = 0,
    RFC_RP_DAMAGE_CALC_METHOD_ELEMENTAR   = 1,
    RFC_RP_DAMAGE_CALC_METHOD_MODIFIED    = 2,
    RFC_RP_DAMAGE_CALC_METHOD_CONSEQUENT  = 3,
};

enum rfc_lc_count_method
{
    RFC_LC_COUNT_METHOD_SLOPES_UP    = 0,
    RFC_LC_COUNT_METHOD_SLOPES_DOWN  = 1,
    RFC_LC_COUNT_METHOD_SLOPES_ALL   = 2,
};
#endif /*!RFC_MINIMAL*/


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
#if RFC_TP_SUPPORT
    RFC_ERROR_TP                    =  4,                           /**< Error while processing turning points */
#endif /*RFC_TP_SUPPORT*/
#if RFC_AT_SUPPORT
    RFC_ERROR_AT                    =  5,                           /**< Error while amplitude transformation */
#endif /*RFC_AT_SUPPORT*/
#if RFC_DH_SUPPORT
    RFC_ERROR_DH_BAD_STREAM         =  6,                           /**< Input stream must be unique */
    RFC_ERROR_DH                    =  7,                           /**< Error while damage history calculation/access */
#endif /*RFC_DH_SUPPORT*/
#if RFC_DAMAGE_FAST
    RFC_ERROR_LUT                   =  8,                           /**< Error while accessing look up tables */
#endif /*RFC_DAMAGE_FAST*/
    RFC_ERROR_DATA_OUT_OF_RANGE     =  9,                           /**< Input data leaves classrange */
    RFC_ERROR_DATA_INCONSISTENT     =  10,                          /**< Processed data is inconsistent (internal error) */
};


#if !RFC_MINIMAL
enum rfc_counting_method
{
    RFC_COUNTING_METHOD_DELEGATED   = -1,                           /**< Method must be implemented via delegator, see member cycle_find_fcn */
    RFC_COUNTING_METHOD_NONE        =  0,                           /**< No counting */
    RFC_COUNTING_METHOD_4PTM        =  1,                           /**< 4 point algorithm (default) */
#if RFC_HCM_SUPPORT
    RFC_COUNTING_METHOD_HCM         =  2,                           /**< 3 point algorithm, Clormann/Seeger (HCM) method */
#endif /*RFC_HCM_SUPPORT*/
#if RFC_ASTM_SUPPORT
    RFC_COUNTING_METHOD_ASTM        =  3,                           /**< 3 point algorithm,  ASTM E1049-85 */
#endif /*RFC_ASTM_SUPPORT*/
    RFC_COUNTING_METHOD_COUNT                                       /**< Number of options */
};
#endif /*!RFC_MINIMAL*/


enum rfc_res_method
{
    /* Don't change order! */
    RFC_RES_NONE                    = 0,                            /**< No residual method */
    RFC_RES_IGNORE                  = 1,                            /**< Ignore residue (same as RFC_RES_NONE) */
    RFC_RES_NO_FINALIZE             = 2,                            /**< Don't finalize the data stream */
#if !RFC_MINIMAL
    RFC_RES_DISCARD                 = 3,                            /**< Discard residue (empty residue) */
    RFC_RES_HALFCYCLES              = 4,                            /**< Related to ASTM, count as half cycles */
    RFC_RES_FULLCYCLES              = 5,                            /**< Count half cycles as full cycles */
    RFC_RES_CLORMANN_SEEGER         = 6,                            /**< Clormann/Seeger method */
    RFC_RES_REPEATED                = 7,                            /**< Repeat residue and count closed cycles */
    RFC_RES_RP_DIN45667             = 8,                            /**< Count residue according to range pair in DIN-45667 */
#endif /*!RFC_MINIMAL*/
    RFC_RES_COUNT                                                   /**< Number of options */
};


#if RFC_DH_SUPPORT
enum rfc_sd_method
{
    RFC_SD_NONE                     = -1,                           /**< No spread damage calculation */
    RFC_SD_HALF_23                  =  0,                           /**< Equally split damage between P2 and P3 */
    RFC_SD_RAMP_AMPLITUDE_23        =  1,                           /**< Spread damage according to amplitude over points between P2 and P3 */
    RFC_SD_RAMP_DAMAGE_23           =  2,                           /**< Spread damage linearly over points between P2 and P3 */
    RFC_SD_RAMP_AMPLITUDE_24        =  3,                           /**< Spread damage according to amplitude over points between P2 and P4 */  
    RFC_SD_RAMP_DAMAGE_24           =  4,                           /**< Spread damage linearly over points between P2 and P4 */
    RFC_SD_FULL_P2                  =  5,                           /**< Assign damage to P2 */
    RFC_SD_FULL_P3                  =  6,                           /**< Assign damage to P3 */
    RFC_SD_TRANSIENT_23             =  7,                           /**< Spread damage transient according to amplitude over points between P2 and P3 */
    RFC_SD_TRANSIENT_23c            =  8,                           /**< Spread damage transient according to amplitude over points between P2 and P4 only until cycle is closed */
    RFC_SD_COUNT                                                    /**< Number of options */
};
#endif /*RFC_DH_SUPPORT*/


enum rfc_wl_defaults
{
    RFC_WL_SD_DEFAULT               =  1000,                        /**< Fatigue strength amplitude (Miner original) */
    RFC_WL_ND_DEFAULT               =  10000000L,                   /**< Cycles according to wl_sd */
    RFC_WL_K_DEFAULT                = -5,                           /**< Woehler slope, always negative */
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
#if !RFC_MINIMAL
typedef     enum        rfc_counting_method     rfc_counting_method_e;      /** Counting method, see RFC_COUNTING... */
typedef     enum        rfc_rp_damage_method    rfc_rp_damage_method_e;     /** Method when calculating damage from range pair counting, see RFC_RP_DAMAGE_CALC_METHOD... */
typedef     enum        rfc_lc_count_method     rfc_lc_count_method_e;      /** Controls which slopes to take into account, when doing the level crossing counting */
#if RFC_DH_SUPPORT
typedef     enum        rfc_sd_method           rfc_sd_method_e;            /** Spread damage method, see RFC_SD... */
#endif /*RFC_DH_SUPPORT*/
typedef     struct      rfc_class_param         rfc_class_param_s;          /** Class parameters (width, offset, count) */
typedef     struct      rfc_wl_param            rfc_wl_param_s;             /** Woehler curve parameters (sd, nd, k, k2, omission) */
typedef     struct      rfc_rfm_item            rfc_rfm_item_s;             /** Rainflow matrix element */
#endif /*!RFC_MINIMAL*/

/* Memory allocation functions typedef */
typedef     void *   ( *rfc_mem_alloc_fcn_t )   ( void *, size_t num, size_t size, int aim );     /** Memory allocation functor */

/* Core functions */
bool        RFC_init                    (       void *ctx, unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                                                           rfc_value_t hysteresis, rfc_flags_e flags );
rfc_state_e RFC_state_get               ( const void *ctx );
rfc_error_e RFC_error_get               ( const void *ctx );
bool        RFC_wl_init_elementary      (       void *ctx, double sx, double nx, double k );
#if !RFC_MINIMAL
bool        RFC_wl_init_original        (       void *ctx, double sd, double nd, double k );
bool        RFC_wl_init_modified        (       void *ctx, double sx, double nx, double k, double k2 );
bool        RFC_wl_init_any             (       void *ctx, const rfc_wl_param_s* wl_param );
bool        RFC_clear_counts            (       void *ctx );
#endif /*!RFC_MINIMAL*/
bool        RFC_deinit                  (       void *ctx );
bool        RFC_feed                    (       void *ctx, const rfc_value_t* data, size_t count );
#if !RFC_MINIMAL
bool        RFC_cycle_process_counts    (       void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags );
bool        RFC_feed_scaled             (       void *ctx, const rfc_value_t* data, size_t count, double factor );
bool        RFC_feed_tuple              (       void *ctx, rfc_value_tuple_s *data, size_t count );
#endif /*!RFC_MINIMAL*/
bool        RFC_finalize                (       void *ctx, rfc_res_method_e residual_method );
#if !RFC_MINIMAL
/* Functions on rainflow matrix */
bool        RFC_rfm_make_symmetric      (       void *ctx );
bool        RFC_rfm_non_zeros           ( const void *ctx, unsigned *count );
bool        RFC_rfm_get                 ( const void *ctx, rfc_rfm_item_s **buffer, unsigned *count );
bool        RFC_rfm_set                 (       void *ctx, const rfc_rfm_item_s *buffer, unsigned count, bool add_only );
bool        RFC_rfm_peek                ( const void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *counts );
bool        RFC_rfm_poke                (       void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t counts, bool add_only );
bool        RFC_rfm_sum                 ( const void *ctx, unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count );
bool        RFC_rfm_damage              ( const void *ctx, unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage );
bool        RFC_rfm_check               ( const void *ctx );
bool        RFC_rfm_refeed              (       void *ctx, rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param );
bool        RFC_lc_get                  ( const void *ctx, rfc_counts_t *lc, rfc_value_t *level );
bool        RFC_lc_from_rfm             ( const void *ctx, rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags );
bool        RFC_lc_from_residue         ( const void *ctx, rfc_counts_t *lc, rfc_value_t *level, rfc_flags_e flags );
bool        RFC_rp_get                  ( const void *ctx, rfc_counts_t *rp, rfc_value_t *Sa );
bool        RFC_rp_from_rfm             ( const void *ctx, rfc_counts_t *rp, rfc_value_t *Sa, const rfc_counts_t *rfm );
bool        RFC_damage                  ( const void *ctx, rfc_value_t *damage, rfc_value_t *damage_residue );
bool        RFC_damage_from_rp          ( const void *ctx, const rfc_counts_t *counts, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_type );
bool        RFC_damage_from_rfm         ( const void *ctx, const rfc_counts_t *rfm, double *damage );
bool        RFC_wl_calc_sx              ( const void *ctx, double s0, double n0, double k, double *sx, double nx, double  k2, double  sd, double nd );
bool        RFC_wl_calc_sd              ( const void *ctx, double s0, double n0, double k, double  sx, double nx, double  k2, double *sd, double nd );
bool        RFC_wl_calc_k2              ( const void *ctx, double s0, double n0, double k, double  sx, double nx, double *k2, double  sd, double nd );
bool        RFC_wl_calc_sa              ( const void *ctx, double s0, double n0, double k, double  n,  double *sa );
bool        RFC_wl_calc_n               ( const void *ctx, double s0, double n0, double k, double  sa, double *n );
bool        RFC_wl_param_set            (       void *ctx, const rfc_wl_param_s * );
bool        RFC_wl_param_get            ( const void *ctx, rfc_wl_param_s * );
bool        RFC_class_param_set         (       void *ctx, const rfc_class_param_s * );
bool        RFC_class_param_get         ( const void *ctx, rfc_class_param_s * );
bool        RFC_class_number            ( const void *ctx, rfc_value_t value, unsigned *class_number );
bool        RFC_class_mean              ( const void *ctx, unsigned class_number, rfc_value_t *class_mean );
bool        RFC_class_upper             ( const void *ctx, unsigned class_number, rfc_value_t *class_upper );
bool        RFC_class_count             ( const void *ctx, unsigned *class_count );
bool        RFC_class_offset            ( const void *ctx, rfc_value_t *class_offset );
bool        RFC_class_width             ( const void *ctx, rfc_value_t *class_width );
bool        RFC_hysteresis              ( const void *ctx, rfc_value_t *hysteresis );
bool        RFC_flags_set               (       void *ctx, int flags, int stack, bool overwrite );
bool        RFC_flags_unset             (       void *ctx, int flags, int stack );
bool        RFC_flags_get               ( const void *ctx, int *flags, int stack );
bool        RFC_flags_check             ( const void *ctx, int flags_to_check, int stack );
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
bool        RFC_tp_init                 (       void *ctx, rfc_value_tuple_s *tp, size_t tp_cap, bool is_static );
bool        RFC_tp_init_autoprune       (       void *ctx, bool autoprune, size_t size, size_t threshold );
bool        RFC_tp_prune                (       void *ctx, size_t count, rfc_flags_e flags );
bool        RFC_tp_refeed               (       void *ctx, rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param );
bool        RFC_tp_clear                (       void *ctx );
#endif /*RFC_TP_SUPPORT*/
bool        RFC_res_get                 ( const void *ctx, const rfc_value_tuple_s **residue, unsigned *count );
#if RFC_DH_SUPPORT
bool        RFC_dh_init                 (       void *ctx, rfc_sd_method_e method, double *dh, size_t dh_cap, bool is_static );
bool        RFC_dh_get                  ( const void *ctx, const double **dh, size_t *count );
#endif /*RFC_DH_SUPPORT*/

#if RFC_AT_SUPPORT
bool        RFC_at_init                 (       void *ctx, const double *Sa, const double *Sm, unsigned count, 
                                                           double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric );
bool        RFC_at_transform            ( const void *ctx, double Sa, double Sm, double *Sa_transformed );
#endif /*RFC_AT_SUPPORT*/

#if RFC_DEBUG_FLAGS
int         RFC_debug_fprintf           (       void *ctx, FILE *stream, const char *fmt, ... );
#endif /*RFC_DEBUG_FLAGS*/


#if RFC_USE_DELEGATES
/* Delegates typedef */
typedef  void                       ( *rfc_cycle_find_fcn_t )    ( rfc_ctx_s *, rfc_flags_e flags );
typedef  bool                       ( *rfc_damage_calc_fcn_t )   ( rfc_ctx_s *, unsigned from_class, unsigned to_class, double *damage, double *Sa_ret );
typedef  bool                       ( *rfc_finalize_fcn_t )      ( rfc_ctx_s *, int residual_methods );
typedef  rfc_value_tuple_s *        ( *rfc_tp_next_fcn_t )       ( rfc_ctx_s *, const rfc_value_tuple_s * );
#if RFC_TP_SUPPORT
typedef  bool                       ( *rfc_tp_set_fcn_t )        ( rfc_ctx_s *, size_t tp_pos, rfc_value_tuple_s * );
typedef  bool                       ( *rfc_tp_get_fcn_t )        ( rfc_ctx_s *, size_t tp_pos, rfc_value_tuple_s ** );
typedef  bool                       ( *rfc_tp_inc_damage_fcn_t ) ( rfc_ctx_s *, size_t tp_pos, double damage );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
typedef  void                       ( *rfc_spread_damage_fcn_t ) ( rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, rfc_flags_e flags );
#endif /*RFC_DH_SUPPORT*/
#if RFC_AT_SUPPORT
typedef  bool                       ( *rfc_at_transform_fcn_t )  ( rfc_ctx_s *, double Sa, double Sm, double *Sa_transformed );
#endif /*RFC_AT_SUPPORT*/
#if RFC_DEBUG_FLAGS
typedef  int                        ( *rfc_debug_vfprintf_fcn_t )( void *, FILE *stream, const char *fmt, va_list arg );
#endif /*RFC_DEBUG_FLAGS*/
#endif /*RFC_USE_DELEGATES*/


/* Value info struct */
struct rfc_value_tuple
{
    rfc_value_t                         value;                      /**< Value. Don't change order, value field must be first! */
    unsigned                            cls;                        /**< Class number, base 0 */
    size_t                              pos;                        /**< Absolute position in input data stream, base 1 */
#if RFC_TP_SUPPORT
    size_t                              adj_pos;                    /**< Absolute position in input data stream of adjacent turning point, base 1. Valid only, if RFC_FLAGS_COUNT_DAMAGE is set! */
    size_t                              tp_pos;                     /**< Position in tp storage, base 1. Only used in residue, in tp storage always 0! */
    rfc_value_t                         avrg;                       /**< Average value of two paired turning points */
#if RFC_DH_SUPPORT
    double                              damage;                     /**< Damage accumulated to this turning point */
#endif /*RFC_DH_SUPPORT*/    
#endif /*RFC_TP_SUPPORT*/
};

#if !RFC_MINIMAL
struct rfc_class_param
{
    unsigned                            count;                      /**< Class count */
    rfc_value_t                         width;                      /**< Class width */
    rfc_value_t                         offset;                     /**< Lower bound of first class */
};

/**
 *   Woehler models
 *  
 *                                \   (S0,N0,k, any point on slope k)
 *                                 \ /
 *                                  o
 *                                   \
 *                                    \
 *                                     \
 *                                      \
 *                                       \   (SX,NX,k2)
 *                                        \ /
 *                                         o
 *                                          \    .             (SD,ND)
 *                                           \         .      /
 *                                            \              o------------------------------ (Miner original)
 *                                             \                   . 
 *                                              \                        . 
 *                                               \                             . 
 *                                                \                                  . 
 *                                                 \                                       . (No fatigue strength) 
 *                                                  (Miner elementar, k2==k)
 *  
 *  ______________________________________________________________________________________ (Omission level)
 *  
 *  s0^|k|  * n0 = sx^|k|  * nx     --->  ln(s0)*|k|  + ln(n0) = ln(sx)*|k|  + ln(nx)
 *  sd^|k2| * nd = sx^|k2| * nx     --->  ln(sd)*|k2| + ln(nd) = ln(sx)*|k2| + ln(nx)
 *  
 *     ln(s0)*|k| + ln(n0) - ln(sd)*|k2|       - ln(nd)                          = ln(sx)*|k| - ln(sx)*|k2|
 *     ln(s0)*|k| + ln(n0) - ln(sd)*|k2|       - ln(nd)                          = ln(sx)*(|k|-|k2|)
 *   ( ln(s0)*|k| + ln(n0) - ln(sd)*|k2|       - ln(nd) ) / (|k|-|k2|)           = ln(sx)
 *   ( ln(s0)*|k| + ln(n0) - ln(sx)*(|k|-|k2|) - ln(nd) ) / |k2|                 = ln(sd)
 *  |( ( ln(s0) - ln(sx) )*|k| + ln(n0)        - ln(nd) ) / ( ln(sd) - ln(sx) )| = |k2|
 */
struct rfc_wl_param
{
    double                              sd;                         /**< Fatigue strength amplitude (Miner original only!) */
    double                              nd;                         /**< Cycles according to wl_sd */
    double                              k;                          /**< Woehler slope, always negative */
    double                              sx;                         /**< Junction point between k and k2 */
    double                              nx;                         /**< Junction point between k and k2 */
    double                              k2;                         /**< Woehler slope between wl_sx and wl_sd, always negative */
    double                              omission;                   /**< Omission level threshold, smaller amplitudes get discarded */
    double                              q;                          /**< Parameter q based on k for "Miner consequent" approach, always positive */
    double                              q2;                         /**< Parameter q based on k2 for "Miner consequent" approach, always positive */
    double                              D;                          /**< If D > 0, parameters define Woehler curve for impaired part */
};

struct rfc_rfm_item
{
    unsigned                            from;                       /**< Start class, base 0 */
    unsigned                            to;                         /**< Ending class, base 0 */
    rfc_counts_t                        counts;                     /**< Counts */
};
#endif /*!RFC_MINIMAL*/


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
#if !RFC_MINIMAL
    rfc_counting_method_e               counting_method;            /**< Searching closed cycles method */
#endif /*!RFC_MINIMAL*/
    rfc_res_method_e                    residual_method;            /**< Used on finalizing */
#if RFC_DH_SUPPORT
    rfc_sd_method_e                     spread_damage_method;       /**< How to spread damage over turning points or time history */
#endif /*RFC_DH_SUPPORT*/

    /* Memory allocation functions */
    rfc_mem_alloc_fcn_t                 mem_alloc;                  /**< Allocate initialized memory */

    /* Counter increments */
    rfc_counts_t                        full_inc;                   /**< Increment for a full cycle */
    rfc_counts_t                        half_inc;                   /**< Increment for a half cycle, used by some residual algorithms */
    rfc_counts_t                        curr_inc;                   /**< Current increment, used by counting algorithms (changed by finalize_res_weight_cycles() only) */

    /* Rainflow class parameters */
    unsigned                            class_count;                /**< Class count */
    rfc_value_t                         class_width;                /**< Class width */
    rfc_value_t                         class_offset;               /**< Lower bound of first class */
    rfc_value_t                         hysteresis;                 /**< Hysteresis filtering, slope must exceed hysteresis to be counted! */

    /* Woehler curve */
    double                              wl_sx;                      /**< Sa of any point on the Woehler curve */
    double                              wl_nx;                      /**< Cycles for Sa on the Woehler curve */
    double                              wl_k;                       /**< Woehler slope, always negative */
#if !RFC_MINIMAL
    double                              wl_sd;                      /**< Fatigue strength amplitude (Miner original) */
    double                              wl_nd;                      /**< Cycles according to wl_sd */
    double                              wl_k2;                      /**< Woehler gradient between wl_sx and wl_sd */
    double                              wl_omission;                /**< Omission level threshold, smaller amplitudes get discarded */
    double                              wl_q;                       /**< Parameter q based on k for "Miner consequent" approach, always positive */
    double                              wl_q2;                      /**< Parameter q based on k2 for "Miner consequent" approach, always positive */
#endif /*!RFC_MINIMAL*/

#if RFC_USE_DELEGATES
    /* Delegates (optional, may be NULL) */
    rfc_tp_next_fcn_t                   tp_next_fcn;                /**< Test for new turning point */
#if RFC_TP_SUPPORT
    rfc_tp_set_fcn_t                    tp_set_fcn;                 /**< Set new turning points */
    rfc_tp_get_fcn_t                    tp_get_fcn;                 /**< Get turning point reference */
    rfc_tp_inc_damage_fcn_t             tp_inc_damage_fcn;          /**< Increase damage for existing turning point */
#endif /*RFC_TP_SUPPORT*/
    rfc_finalize_fcn_t                  finalize_fcn;               /**< Finalizing function */
    rfc_cycle_find_fcn_t                cycle_find_fcn;             /**< Find next cycle(s) and process */
    rfc_damage_calc_fcn_t               damage_calc_fcn;            /**< Damage calculating function */
#if RFC_DH_SUPPORT
    rfc_spread_damage_fcn_t             spread_damage_fcn;          /**< Spread damage over turning points and damage history */
#endif /*RFC_DH_SUPPORT*/
#if RFC_AT_SUPPORT
    rfc_at_transform_fcn_t              at_transform_fcn;           /**< Amplitude transformation to take mean load influence into account */
#endif /*RFC_AT_SUPPORT*/
#if RFC_DEBUG_FLAGS
    rfc_debug_vfprintf_fcn_t            debug_vfprintf_fcn;         /**< Debug vfprintf() function */
#endif /*RFC_DEBUG_FLAGS*/
#endif /*RFC_USE_DELEGATES*/
    
    /* Residue */
    rfc_value_tuple_s                  *residue;                    /**< Buffer for residue */
    size_t                              residue_cap;                /**< Buffer capacity in number of elements (max. 2*class_count) */
    size_t                              residue_cnt;                /**< Number of value tuples in buffer */

    /* Non-sparse storages (optional, may be NULL) */
    rfc_counts_t                       *rfm;                        /**< Rainflow matrix, always class_count^2 elements (row-major, row=from, to=col). */
#if !RFC_MINIMAL
    rfc_counts_t                       *rp;                         /**< Range pair counts, always class_count elements */
    rfc_counts_t                       *lc;                         /**< Level crossing counts, always class_count elements. Every per .flags selected slope increments by .full_inc! */
#endif /*!RFC_MINIMAL*/

#if RFC_TP_SUPPORT
    /* Turning points storage (optional, may be NULL) */
    rfc_value_tuple_s                  *tp;                         /**< Buffer for turning points, pointer may be changed whilst memory reallocation! */
    size_t                              tp_cap;                     /**< Buffer capacity (number of elements) */
    size_t                              tp_cnt;                     /**< Number of turning points in buffer */
    int                                 tp_locked;                  /**< If tp_locked > 0, no more points can be added. RFC_tp_prune() may delete content. Field .damage is always mutable */
    size_t                              tp_prune_size;              /**< Size for autoprune */
    size_t                              tp_prune_threshold;         /**< Threshold for (auto)pruning */
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    const rfc_value_t                  *dh_istream;                 /**< Input stream */
    double                             *dh;                         /**< Damage history, pointer may be changed whilst memory reallocation! */
    size_t                              dh_cap;                     /**< Capacity of dh */
    size_t                              dh_cnt;                     /**< Number of values in dh */
#endif /*RFC_DH_SUPPORT*/

    /* Damage */
#if RFC_DAMAGE_FAST
    double                             *damage_lut;                 /**< Damage look-up table */
    int                                 damage_lut_inapt;           /**< Greater 0, if values in damage_lut aren't proper to Woehler curve parameters */
#if RFC_AT_SUPPORT
    double                             *amplitude_lut;              /**< Amplitude look-up table, only valid if damage_lut_inapt == 0 */
#endif /*RFC_AT_SUPPORT*/
#endif /*RFC_DAMAGE_FAST*/
    double                              damage;                     /**< Cumulated damage (damage resulting from residue included) */
    double                              damage_residue;             /**< Partial damage in .damage influenced by taking residue into account (after finalizing) */

#if RFC_AT_SUPPORT
    struct at
    {
        const double                   *Sa;
        const double                   *Sm;
        unsigned                        count;
        double                          M;
        double                          Sm_rig;
        double                          R_rig;
        bool                            R_pinned;
    }                                   
                                        at;
#endif /*RFC_AT_SUPPORT*/

    /* Internal usage */
    struct internal
    {
        int                             flags;                      /**< Flags (enum rfc_flags) */
#if _DEBUG
        bool                            finalizing;                 /**< true, when finalizing */
#endif /*_DEBUG*/
#if RFC_DEBUG_FLAGS
        int                             debug_flags;                /**< Flags for debugging */
#endif /*RFC_DEBUG_FLAGS*/
        int                             slope;                      /**< Current signal slope */
        rfc_value_tuple_s               extrema[2];                 /**< Local or global extrema depending on RFC_GLOBAL_EXTREMA */
#if RFC_GLOBAL_EXTREMA
        bool                            extrema_changed;            /**< True if one extrema has changed */
#endif /*!RFC_GLOBAL_EXTREMA*/
        size_t                          pos;                        /**< Absolute position in data input stream, base 1 */
        size_t                          pos_offset;                 /**< Offset for pos */
        rfc_value_tuple_s               residue[3];                 /**< Static residue (if class_count is zero) */
        size_t                          residue_cap;                /**< Capacity of static residue */
        bool                            res_static;                 /**< true, if .residue refers the static residue .internal.residue */
#if !RFC_MINIMAL
        rfc_wl_param_s                  wl;                         /**< Shadowed Woehler curve parameters */
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
        rfc_value_tuple_s               margin[2];                  /**< First and last data point */
        int                             margin_stage;               /**< 0: Init, 1: Left margin set, 2: 1st turning point is safe */
        bool                            tp_static;                  /**< true, if tp is statically allocated */
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
        bool                            dh_static;                  /**< true, if dh is statically allocated */
#endif /*RFC_DH_SUPPORT*/
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
#if RFC_AT_SUPPORT
        struct
        {
            double                      Sa[5];
            double                      Sm[5];
            unsigned                    count;
        }                               at_haigh;
#endif /*RFC_AT_SUPPORT*/
#if RFC_USE_DELEGATES
            void                       *obj;
#endif /*RFC_USE_DELEGATES*/
    }
                                        internal;
};

#ifdef __cplusplus
}  /* extern "C" */
}  /* namespace RFC_CPP_NAMESPACE */
#endif /*__cplusplus*/

#pragma pack(pop)

#endif /*RAINFLOW_H*/
