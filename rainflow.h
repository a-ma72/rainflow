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
 * []  "Schaedigungsbasierte Hysteresefilter"; Hack, M, D386 (Diss Univ. Kaiserslautern), Shaker Verlag Aachen, 1998, ISBN 3-8265-3936-2
 * []  "Hysteresis and Phase Transition"
 *     Brokate, M.; Sprekels, J.; Applied Mathematical Sciences 121; Springer, New York, 1996
 * []  "Rainflow counting and energy dissipation in elastoplasticity"; Eur. J. Mech. A/Solids 15, pp. 705-737, 1996
 *     Brokate, M.; Dressler, K.; Krejci, P.
 * []  "Multivariate Density Estimation: Theory, Practice and Visualization". New York, Chichester, Wiley & Sons, 1992
 *     Scott, D.
 * []  "Werkstoffmechanik - Bauteile sicher beurteilen undWerkstoffe richtig einsetzen"; 
 *      Ralf Bürgel, Hans Albert Richard, Andre Riemer; Springer FachmedienWiesbaden 2005, 2014
 * [] "Zählverfahren und Lastannahme in der Betriebsfestigkeit";
 *    Michael Köhler, Sven Jenne • Kurt Pötter, Harald Zenner; Springer-Verlag Berlin Heidelberg 2012
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

#if RFC_MINIMAL
#undef  RFC_USE_DELEGATES
#define RFC_USE_DELEGATES    OFF
#undef  RFC_HCM_SUPPORT
#define RFC_HCM_SUPPORT      OFF
#undef  RFC_TP_SUPPORT
#define RFC_TP_SUPPORT       OFF
#undef  RFC_DH_SUPPORT
#define RFC_DH_SUPPORT       OFF
#undef  RFC_AT_SUPPORT
#define RFC_AT_SUPPORT       OFF
#undef  RFC_GLOBAL_EXTREMA
#define RFC_GLOBAL_EXTREMA   OFF
#undef  RFC_DAMAGE_FAST
#define RFC_DAMAGE_FAST      OFF
#else /*!RFC_MINIMAL*/
#ifndef RFC_USE_DELEGATES
#define RFC_USE_DELEGATES OFF
#endif /*RFC_USE_DELEGATES*/
#ifndef RFC_HCM_SUPPORT
#define RFC_HCM_SUPPORT OFF
#endif /*RFC_HCM_SUPPORT*/
#ifndef RFC_TP_SUPPORT
#define RFC_TP_SUPPORT ON
#endif /*RFC_TP_SUPPORT*/
#ifndef RFC_DH_SUPPORT
#define RFC_DH_SUPPORT OFF
#endif /*RFC_DH_SUPPORT*/
#ifndef RFC_AT_SUPPORT
#define RFC_AT_SUPPORT OFF
#endif /*RFC_AT_SUPPORT*/
#ifndef RFC_GLOBAL_EXTREMA
#define RFC_GLOBAL_EXTREMA OFF
#endif /*RFC_GLOBAL_EXTREMA*/
#ifndef RFC_DAMAGE_FAST
#define RFC_DAMAGE_FAST ON
#endif /*RFC_DAMAGE_FAST*/
#endif /*RFC_MINIMAL*/

/* Memory allocation aim info */
enum
{
    RFC_MEM_AIM_TEMP                =  0,
    RFC_MEM_AIM_RESIDUE             =  1,
    RFC_MEM_AIM_MATRIX              =  2,
    RFC_MEM_AIM_RP                  =  3,
    RFC_MEM_AIM_LC                  =  4,
#if RFC_TP_SUPPORT
    RFC_MEM_AIM_TP                  =  5,
#endif
#if RFC_DAMAGE_FAST
    RFC_MEM_AIM_DLUT                =  6,
#endif
#if RFC_HCM_SUPPORT
    RFC_MEM_AIM_HCM                 =  7,
#endif
#if RFC_DH_SUPPORT
    RFC_MEM_AIM_DH                  =  8,
#endif
};


/* Memory allocation functions typedef */
typedef void * ( *rfc_mem_alloc_fcn_t )( void *, size_t num, size_t size, int aim );

/* Typedefs */
typedef RFC_VALUE_TYPE          RFC_value_type;      /** Input data value type */
typedef RFC_COUNTS_VALUE_TYPE   RFC_counts_type;     /** Type of counting values */
typedef struct rfc_ctx          rfc_ctx_s;           /** Forward declaration (rainflow context) */
typedef struct rfc_value_tuple  rfc_value_tuple_s;   /** Tuple of value and index position */
#if !RFC_MINIMAL
typedef struct rfc_class_param  rfc_class_param_s;   /** Class parameters (width, offset, count) */
#endif /*RFC_MINIMAL*/


/* Core functions */
bool   RFC_init              ( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                                          RFC_value_type hysteresis );
bool   RFC_deinit            ( void *ctx );
int    RFC_flags_set         ( void *ctx, int flags );
bool   RFC_feed              ( void *ctx, const RFC_value_type* data, size_t count );
#if !RFC_MINIMAL
bool   RFC_feed_tuple        ( void *ctx, rfc_value_tuple_s *data, size_t count );
#endif /*RFC_MINIMAL*/
bool   RFC_finalize          ( void *ctx, int residual_method );

#if RFC_TP_SUPPORT
bool   RFC_tp_init           ( void *ctx, rfc_value_tuple_s *tp, size_t tp_cap, bool is_static );
bool   RFC_tp_init_autoprune ( void *ctx, bool autoprune, size_t size, size_t threshold );
bool   RFC_tp_prune          ( void *ctx, size_t count, int flags );
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
bool   RFC_dh_init           ( void *ctx, double *dh, size_t dh_cap, bool is_static );
#endif /*RFC_DH_SUPPORT*/

#if RFC_AT_SUPPORT
bool   RFC_at_init           ( void *ctx, const double *Sa, const double *Sm, unsigned count, 
                                          double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric );
double RFC_at_transform      ( void *ctx, double Sa, double Sm );
#endif /*RFC_AT_SUPPORT*/

#if RFC_USE_DELEGATES
/* Delegates typedef */
typedef  void                ( *rfc_cycle_find_fcn_t )    ( rfc_ctx_s *, int flags );
typedef  double              ( *rfc_damage_calc_fcn_t )   ( rfc_ctx_s *, unsigned from_class, unsigned to_class );
typedef  bool                ( *rfc_finalize_fcn_t )      ( rfc_ctx_s *, int residual_methods );
typedef  rfc_value_tuple_s * ( *rfc_tp_next_fcn_t )       ( rfc_ctx_s *, const rfc_value_tuple_s * );
#if RFC_TP_SUPPORT
typedef  bool                ( *rfc_tp_add_fcn_t )        ( rfc_ctx_s *, rfc_value_tuple_s * );
typedef  bool                ( *rfc_tp_prune_fcn_t )      ( rfc_ctx_s *, size_t, int );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
typedef  void                ( *rfc_spread_damage_fcn_t ) ( rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, int flags );
#endif /*RFC_DH_SUPPORT*/
#if RFC_AT_SUPPORT
typedef  double              ( *rfc_at_transform_fcn_t )  ( rfc_ctx_s *, double Sa, double Sm );
#endif /*RFC_AT_SUPPORT*/
#endif /*RFC_USE_DELEGATES*/


/* Value info struct */
typedef struct rfc_value_tuple
{
    RFC_value_type                      value;                      /**< Value. Don't change order, value field must be first! */
    unsigned                            class;                      /**< Class number, base 0 */
    size_t                              pos;                        /**< Absolute position in input data stream, base 1 */
#if RFC_TP_SUPPORT
    size_t                              tp_pos;                     /**< Position in tp storage, base 1 */
#if RFC_DH_SUPPORT
    double                              damage;                     /**< Damage accumulated to this turning point */
#endif /*RFC_DH_SUPPORT*/    
#endif /*RFC_TP_SUPPORT*/
} rfc_value_tuple_s;

#if !RFC_MINIMAL
typedef struct rfc_class_param
{
    unsigned                            count;                      /**< Class count */
    RFC_value_type                      width;                      /**< Class width */
    RFC_value_type                      offset;                     /**< Lower bound of first class */
} rfc_class_param_s;
#endif /*RFC_MINIMAL*/


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
        RFC_ERROR_NOERROR,                                          /**< No error */
        RFC_ERROR_INVARG,                                           /**< Invalid arguments passed */
        RFC_ERROR_MEMORY,                                           /**< Error on memory allocation */
    }                                   error;                      /**< Error code */

    enum
    {
        RFC_FLAGS_COUNT_MATRIX          = 1 << 0,                   /**< Count into matrix */
        RFC_FLAGS_COUNT_DAMAGE          = 1 << 1,                   /**< Count pseudo damage */
#if !RFC_MINIMAL
#if RFC_DH_SUPPORT
        RFC_FLAGS_COUNT_DH              = 1 << 2,                   /**< Spread damage */
#endif /*RFC_DH_SUPPORT*/
        RFC_FLAGS_COUNT_RP              = 1 << 3,                   /**< Count into range pair */
        RFC_FLAGS_COUNT_LC_UP           = 1 << 4,                   /**< Count into level crossing (only rising slopes) */
        RFC_FLAGS_COUNT_LC_DN           = 1 << 5,                   /**< Count into level crossing (only falling slopes) */
        RFC_FLAGS_COUNT_LC              = RFC_FLAGS_COUNT_LC_UP     /**< Count into level crossing (all slopes) */
                                        | RFC_FLAGS_COUNT_LC_DN,
        RFC_FLAGS_ENFORCE_MARGIN        = 1 << 6,                   /**< Enforce first and last data point are turning points */
#endif /*RFC_MINIMAL*/
        RFC_FLAGS_COUNT_ALL             = RFC_FLAGS_COUNT_MATRIX    /**< Count all */
                                        | RFC_FLAGS_COUNT_DAMAGE
#if RFC_DH_SUPPORT
                                        | RFC_FLAGS_COUNT_DH
#endif /*RFC_DH_SUPPORT*/
#if !RFC_MINIMAL
                                        | RFC_FLAGS_COUNT_RP
                                        | RFC_FLAGS_COUNT_LC,
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
        RFC_FLAGS_TPPRUNE_PRESERVE_POS  = 1 << 7,
        RFC_FLAGS_TPPRUNE_PRESERVE_RES  = 1 << 8,
        RFC_FLAGS_TPAUTOPRUNE           = 1 << 9,                   /**< Automatic prune on tp */
#endif /*RFC_TP_SUPPORT*/
    }
                                        flags;                      /**< Flags */
#if !RFC_MINIMAL
    enum
    {
        RFC_COUNTING_METHOD_DELEGATED   = -1,                       /**< Method must be implemented via delegator, see member cycle_find_fcn */
        RFC_COUNTING_METHOD_NONE        =  0,                       /**< No counting */
        RFC_COUNTING_METHOD_4PTM        =  1,                       /**< 4 point algorithm (default) */
#if RFC_HCM_SUPPORT
        RFC_COUNTING_METHOD_HCM         =  2,                       /**< 3 point algorithm, Clormann/Seeger (HCM) method */
#endif /*RFC_HCM_SUPPORT*/
    }
                                        counting_method;
#endif /*RFC_MINIMAL*/

    enum 
    {
        /* Don't change order! */
        RFC_RES_NONE                    = 0,                        /**< No residual method */
        RFC_RES_IGNORE,                                             /**< Ignore residue (same as RFC_RES_NONE) */
#if !RFC_MINIMAL
        RFC_RES_DISCARD,                                            /**< Discard residue (empty residue) */
        RFC_RES_HALFCYCLES,                                         /**< ASTM */
        RFC_RES_FULLCYCLES,                                         /**< Count half cycles as full cycles */
        RFC_RES_CLORMANN_SEEGER,                                    /**< Clormann/Seeger method */
        RFC_RES_REPEATED,                                           /**< Repeat residue and count closed cycles */
        RFC_RES_RP_DIN45667,                                        /**< Count residue according to range pair in DIN-45667 */
#endif /*RFC_MINIMAL*/
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
#if !RFC_MINIMAL
    double                              wl_k2;                      /**< Woehler gradient below wl_sd */
    double                              wl_omission;                /**< Omission level */
#endif /*RFC_MINIMAL*/

#if RFC_USE_DELEGATES
    /* Delegates (optional, may be NULL) */
    rfc_tp_next_fcn_t                   tp_next_fcn;                /**< Test if new turning point exists */
#if RFC_TP_SUPPORT
    rfc_tp_add_fcn_t                    tp_add_fcn;                 /**< Handling new turning points */
    rfc_tp_prune_fcn_t                  tp_prune_fcn;               /**< Prune turning points */
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
#endif /*RFC_USE_DELEGATES*/
    
    /* Residue */
    rfc_value_tuple_s                  *residue;                    /**< Buffer for residue */
    size_t                              residue_cap;                /**< Buffer capacity in number of elements (max. 2*class_count) */
    size_t                              residue_cnt;                /**< Number of value tuples in buffer */

    /* Non-sparse storages (optional, may be NULL) */
    RFC_counts_type                    *matrix;                     /**< Rainflow matrix, always class_count^2 elements */
#if !RFC_MINIMAL
    RFC_counts_type                    *rp;                         /**< Range pair counts, always class_count elements */
    RFC_counts_type                    *lc;                         /**< Level crossing counts, always class_count elements */
#endif /*RFC_MINIMAL*/

#if RFC_TP_SUPPORT
    /* Turning points storage (optional, may be NULL) */
    rfc_value_tuple_s                  *tp;                         /**< Buffer for turning points, pointer may be changed thru memory reallocation! */
    size_t                              tp_cap;                     /**< Buffer capacity (number of elements) */
    size_t                              tp_cnt;                     /**< Number of turning points in buffer */
    bool                                tp_locked;                  /**< If tp_locked, tp is freezed. Only RFC_tp_prune() may change content */
    size_t                              tp_prune_size;              /**< Size for autoprune */
    size_t                              tp_prune_threshold;         /**< Threshold for (auto)pruning */
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    double                             *dh;                         /**< Damage history, pointer may be changed thru memory reallocation! */
    size_t                              dh_cap;                     /**< Capacity of dh */
    size_t                              dh_cnt;                     /**< Number of values in dh */
#endif /*RFC_DH_SUPPORT*/

    /* Damage */
#if RFC_DAMAGE_FAST
    double                             *damage_lut;                 /**< Damage look-up table */
#endif /*RFC_DAMAGE_FAST*/
    double                              pseudo_damage;              /**< Cumulated pseudo damage */

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
        int                             slope;                      /**< Current signal slope */
        rfc_value_tuple_s               extrema[2];                 /**< Local or global extrema depending on RFC_GLOBAL_EXTREMA */
#if !RFC_MINIMAL
        bool                            extrema_changed;            /**< True if one extrema has changed */
#endif /*RFC_MINIMAL*/
        size_t                          pos;                        /**< Absolute position in data input stream, base 1 */
        size_t                          global_offset;              /**< Offset for pos */
        rfc_value_tuple_s               residue[3];                 /**< Static residue (if class_count is zero) */
        size_t                          residue_cap;                /**< Capacity of static residue */
        bool                            res_static;                 /**< true, if static residue is in use */
#if RFC_TP_SUPPORT
        rfc_value_tuple_s               margin[2];                  /**< First and last data point */
        int                             margin_stage;               /**< 0: Init, 1: Left margin set, 2: 1st turning point safe */
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
        }                               at;
#endif /*RFC_AT_SUPPORT*/

    }
                                        internal;
} rfc_ctx_s;
