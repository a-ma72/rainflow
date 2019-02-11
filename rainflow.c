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

#if COAN_INVOKED
/* This version is generated via coan (http://coan2.sourceforge.net/) */
#endif /*COAN_INVOKED*/


#include "rainflow.h"

#include <assert.h>  /* assert() */
#include <math.h>    /* exp(), log(), fabs() */
#include <stdlib.h>  /* calloc(), free(), abs() */
#include <string.h>  /* memset() */
#include <float.h>   /* DBL_MAX */


#if MATLAB_MEX_FILE
#if !RFC_MINIMAL
#define RFC_MEX_USAGE \
"\nUsage:\n"\
"[pd,re,rm,rp,lc,tp] = rfc( 'rfc', data, class_count, class_width, class_offset, hysteresis, residual_method, enfore_margin, use_hcm )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"\
"    rp = Range pair counts\n"\
"    lc = Level crossings\n"\
"    tp = Turning points\n"\
"\n"\
"[Sa] = rfc( 'amptransform', Sa, Sm, M, target, R_pinned )\n"\
"             Sa = Amplitude\n"\
"             Sm = Mean load\n"\
"              M = Mean load sensitivity\n"\
"         target = Mean load or mean load ratio (R)\n"\
"    target_is_R = true, if target is R (otherwise target is test rig mean load)\n"
#else /*RFC_MINIMAL*/
#define RFC_MEX_USAGE \
"\nUsage:\n"\
"[pd,re,rm] = rfc( data, class_count, class_width, class_offset, hysteresis )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"
#endif /*!RFC_MINIMAL*/
#pragma message(RFC_MEX_USAGE)
#include <ctype.h>
#include <string.h>
#include <mex.h>
#endif /*MATLAB_MEX_FILE*/

/* Core functions */
#if !RFC_MINIMAL
#if RFC_AT_SUPPORT
static void                 clear_at                        (       rfc_ctx_s * );
#endif /*RFC_AT_SUPPORT*/
#if RFC_DAMAGE_FAST
static void                 clear_lut                       (       rfc_ctx_s * );
#endif /*RFC_DAMAGE_FAST*/
static void                 cycle_find                      (       rfc_ctx_s *, rfc_flags_e flags );
#else /*RFC_MINIMAL*/
#define cycle_find          cycle_find_4ptm
#endif /*!RFC_MINIMAL*/
static bool                 feed_once                       (       rfc_ctx_s *, const rfc_value_tuple_s* tp, rfc_flags_e flags );
#if RFC_DH_SUPPORT
static bool                 feed_once_dh                    (       rfc_ctx_s *, const rfc_value_tuple_s* pt );
#endif /*RFC_DH_SUPPORT*/
#if RFC_TP_SUPPORT
static bool                 feed_once_tp_check_margin       (       rfc_ctx_s *, const rfc_value_tuple_s* pt, rfc_value_tuple_s** tp_residue );
#endif /*RFC_TP_SUPPORT*/
static bool                 feed_finalize                   (       rfc_ctx_s * );
#if RFC_TP_SUPPORT
static bool                 feed_finalize_tp                (       rfc_ctx_s *, rfc_value_tuple_s *tp_interim, rfc_flags_e flags );
#endif /*RFC_TP_SUPPORT*/
#if RFC_HCM_SUPPORT
static bool                 feed_finalize_hcm               (       rfc_ctx_s *, rfc_flags_e flags );
#endif /*!RFC_HCM_SUPPORT*/
static rfc_value_tuple_s *  feed_filter_pt                  (       rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 cycle_find_4ptm                 (       rfc_ctx_s *, rfc_flags_e flags );
#if RFC_HCM_SUPPORT
static void                 cycle_find_hcm                  (       rfc_ctx_s *, rfc_flags_e flags );
#endif /*RFC_HCM_SUPPORT*/
#if !RFC_MINIMAL
static void                 cycle_process_lc                (       rfc_ctx_s *, rfc_flags_e flags );
#endif /*!RFC_MINIMAL*/
static void                 cycle_process_counts            (       rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, rfc_flags_e flags );
/* Methods on residue */
static bool                 finalize_res_ignore             (       rfc_ctx_s *, rfc_flags_e flags );
#if !RFC_MINIMAL
static bool                 finalize_res_discard            (       rfc_ctx_s *, rfc_flags_e flags );
static bool                 finalize_res_weight_cycles      (       rfc_ctx_s *, rfc_counts_t weight, rfc_flags_e flags );
static bool                 finalize_res_clormann_seeger    (       rfc_ctx_s *, rfc_flags_e flags );
static bool                 RFC_finalize_res_rp_DIN45667    (       rfc_ctx_s *, rfc_flags_e flags );
static bool                 finalize_res_repeated           (       rfc_ctx_s *, rfc_flags_e flags );
static bool                 residue_exchange                (       rfc_ctx_s *, rfc_value_tuple_s **residue, size_t *residue_cap, size_t *residue_cnt, bool restore );
#endif /*!RFC_MINIMAL*/
static void                 residue_remove_item             (       rfc_ctx_s *, size_t index, size_t count );
/* Memory allocator */
static void *               mem_alloc                       ( void *ptr, size_t num, size_t size, rfc_mem_aim_e aim );
#if RFC_TP_SUPPORT
/* Methods on turning points history */
static bool                 tp_set                          (       rfc_ctx_s *, size_t tp_pos, rfc_value_tuple_s *pt );
static bool                 tp_get                          (       rfc_ctx_s *, size_t tp_pos, rfc_value_tuple_s **pt );
static bool                 tp_inc_damage                   (       rfc_ctx_s *, size_t tp_pos, double damage );
static void                 tp_lock                         (       rfc_ctx_s *, bool do_lock );
static bool                 tp_refeed                       (       rfc_ctx_s *, rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
static bool                 spread_damage                   (       rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, rfc_flags_e flags );
#endif /*RFC_DH_SUPPORT*/
#if RFC_AT_SUPPORT
static bool                 at_R_to_Sm_norm                 (       rfc_ctx_s *, double R, double *Sm_norm );
static bool                 at_alleviation                  (       rfc_ctx_s *, double Sm_norm, double *alleviation );
#endif /*RFC_AT_SUPPORT*/
/* Other */
static bool                 damage_calc_amplitude           (       rfc_ctx_s *, double Sa, double *damage );
static bool                 damage_calc                     (       rfc_ctx_s *, unsigned class_from, unsigned class_to, double *damage, double *Sa_ret );
#if RFC_DAMAGE_FAST
static void                 damage_lut_init                 (       rfc_ctx_s * );
static bool                 damage_calc_fast                (       rfc_ctx_s *, unsigned class_from, unsigned class_to, double *damage, double *Sa_ret );
#endif /*RFC_DAMAGE_FAST*/
static bool                 error_raise                     (       rfc_ctx_s *, rfc_error_e );
static rfc_value_t          value_delta                     (       rfc_value_t from_val, rfc_value_t to, int *sign_ptr );


#define QUANTIZE( r, v )    ( (unsigned)( ((v) - (r)->class_offset) / (r)->class_width ) )
#define AMPLITUDE( r, i )   ( (double)(r)->class_width * (i) / 2 )
#define CLASS_MEAN( r, c )  ( (double)(r)->class_width * (0.5 + (c)) + (r)->class_offset )
#define CLASS_UPPER( r, c ) ( (double)(r)->class_width * (1.0 + (c)) + (r)->class_offset )
#define NUMEL( x )          ( sizeof(x) / sizeof(*(x)) )

#define RFC_CTX_CHECK_AND_ASSIGN                                                    \
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;                                           \
                                                                                    \
    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )                         \
    {                                                                               \
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );                            \
    }                                                                               \



#if !RFC_TP_SUPPORT
/**
 * @brief      Initialization (rainflow context).
 *
 * @param      ctx           The rainflow context
 * @param      class_count   The class count
 * @param      class_width   The class width
 * @param      class_offset  The class offset
 * @param      hysteresis    The hysteresis 
 * @param      flags         The flags
 *
 * @return     true on success
 */
#else /*RFC_TP_SUPPORT*/
/**
 * @brief      Initialization (rainflow context).
 *
 * @param      ctx           The rainflow context
 * @param      class_count   The class count
 * @param      class_width   The class width
 * @param      class_offset  The class offset
 * @param      hysteresis    The hysteresis 
 * @param      flags         The flags
 * @param      tp            Pointer to turning points buffer
 * @param      tp_cap        Number of turning points in buffer 
 *
 * @return     true on success
 */
#endif /*!RFC_TP_SUPPORT*/
bool RFC_init( void *ctx, unsigned class_count, rfc_value_t class_width, rfc_value_t class_offset, 
                          rfc_value_t hysteresis, rfc_flags_e flags )
{
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state != RFC_STATE_INIT0 )
    {
        return false;
    }

    /* Flags */
    if( flags == RFC_FLAGS_DEFAULT )
    {
#if !RFC_MINIMAL
        flags                               = RFC_FLAGS_COUNT_ALL            | 
#if _DEBUG
                                              RFC_FLAGS_LOG_CLOSED_CYCLES    |
#endif /*_DEBUG*/
#if RFC_TP_SUPPORT  
                                              RFC_FLAGS_TPPRUNE_PRESERVE_POS | 
                                              RFC_FLAGS_TPPRUNE_PRESERVE_RES |
#endif /*RFC_TP_SUPPORT*/
                                              0;
#else /*RFC_MINIMAL*/
        flags                               = RFC_FLAGS_COUNT_RFM            | 
                                              RFC_FLAGS_COUNT_DAMAGE;
#endif /*!RFC_MINIMAL*/
    }
    rfc_ctx->internal.flags                 = flags;

    /* Counter increments */
    rfc_ctx->full_inc                       = RFC_FULL_CYCLE_INCREMENT;
    rfc_ctx->half_inc                       = RFC_HALF_CYCLE_INCREMENT;
    rfc_ctx->curr_inc                       = RFC_FULL_CYCLE_INCREMENT;

    if( class_count )
    {
        if( class_count > 512 || class_width <= 0.0 )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
    }

    /* Rainflow class parameters */
    rfc_ctx->class_count                    = class_count;
    rfc_ctx->class_width                    = class_width;
    rfc_ctx->class_offset                   = class_offset;
    rfc_ctx->hysteresis                     = hysteresis;

    /* Values for a "pseudo Woehler curve" */
    rfc_ctx->state = RFC_STATE_INIT;   /* Bypass sanity check for state in wl_init() */
    RFC_wl_init_elementary( rfc_ctx, 1e3 /*sx*/, 1e7 /*nx*/, -5.0 /*k*/ );
    rfc_ctx->state = RFC_STATE_INIT0;  /* Reset state */

    /* Memory allocator */
    if( !rfc_ctx->mem_alloc )
    {
        rfc_ctx->mem_alloc = mem_alloc;
    }
    
#if RFC_USE_DELEGATES
    /* Delegates (optional, set to NULL for standard or to your own functions! ) */
#if RFC_TP_SUPPORT
    rfc_ctx->tp_next_fcn                    = NULL;
    rfc_ctx->tp_set_fcn                     = NULL;
    rfc_ctx->tp_get_fcn                     = NULL;
    rfc_ctx->tp_inc_damage_fcn              = NULL;
#endif /*RFC_TP_SUPPORT*/
    rfc_ctx->cycle_find_fcn                 = NULL;
    rfc_ctx->finalize_fcn                   = NULL;
    rfc_ctx->damage_calc_fcn                = NULL;
#endif /*RFC_USE_DELEGATES*/

#if !RFC_MINIMAL
    /* Rainflow counting method */
    rfc_ctx->counting_method                = RFC_COUNTING_METHOD_4PTM;
#endif /*!RFC_MINIMAL*/

    /* Residue */
    rfc_ctx->internal.residue_cap           = NUMEL( rfc_ctx->internal.residue );
    rfc_ctx->residue_cnt                    = 0;
    rfc_ctx->residue_cap                    = 2 * rfc_ctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */

    if( rfc_ctx->residue_cap <= rfc_ctx->internal.residue_cap )
    {
        rfc_ctx->residue_cap                = rfc_ctx->internal.residue_cap; /* At least 3 elements are needed (two to define a slope and one as interim point) */
        rfc_ctx->residue                    = rfc_ctx->internal.residue;
        rfc_ctx->internal.res_static        = true;
    }
    else
    {
        rfc_ctx->residue                    = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, rfc_ctx->residue_cap, 
                                                                                      sizeof(rfc_value_tuple_s), RFC_MEM_AIM_RESIDUE );
        rfc_ctx->internal.res_static        = false;
    }

    if( rfc_ctx->class_count )
    {
        int ok = rfc_ctx->residue != NULL;

        if( ok && ( flags & RFC_FLAGS_COUNT_RFM ) )
        {
            /* Non-sparse storages (optional, may be NULL) */
            rfc_ctx->rfm                    = (rfc_counts_t*)rfc_ctx->mem_alloc( NULL, class_count * class_count, 
                                                                                 sizeof(rfc_counts_t), RFC_MEM_AIM_MATRIX );
            if( !rfc_ctx->rfm ) ok = false;
        }
#if !RFC_MINIMAL
        if( ok && ( flags & RFC_FLAGS_COUNT_RP ) )
        {
            rfc_ctx->rp                     = (rfc_counts_t*)rfc_ctx->mem_alloc( NULL, class_count,
                                                                                 sizeof(rfc_counts_t), RFC_MEM_AIM_RP );
            if( !rfc_ctx->rp ) ok = false;
        }

        if( ok && ( flags & RFC_FLAGS_COUNT_LC ) )
        {
            rfc_ctx->lc                     = (rfc_counts_t*)rfc_ctx->mem_alloc( NULL, class_count,
                                                                                 sizeof(rfc_counts_t), RFC_MEM_AIM_LC );
            if( !rfc_ctx->lc ) ok = false;
        }
#endif /*!RFC_MINIMAL*/
        if( !ok )
        {
            RFC_deinit( rfc_ctx );
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
    }

    /* Damage */
    rfc_ctx->damage                         = 0.0;
    rfc_ctx->damage_residue                 = 0.0;

    /* Internals */
    rfc_ctx->internal.slope                 = 0;
    rfc_ctx->internal.extrema[0]            = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]            = nil;  /* local maximum */
#if RFC_GLOBAL_EXTREMA    
    rfc_ctx->internal.extrema_changed       = false;
#endif /*RFC_GLOBAL_EXTREMA*/
#if !RFC_MINIMAL
    /* Make a shadow copy of the Woehler curve parameters */
    rfc_ctx->state = RFC_STATE_INIT;   /* Bypass sanity check for state in wl_init() */
    RFC_wl_param_get( rfc_ctx, &rfc_ctx->internal.wl );
    rfc_ctx->state = RFC_STATE_INIT0;  /* Reset state */
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
    rfc_ctx->internal.margin[0]             = nil;  /* left  margin */
    rfc_ctx->internal.margin[1]             = nil;  /* right margin */
    rfc_ctx->internal.margin_stage          = 0;
    /* Turning points storage (optional, may be NULL) */
    rfc_ctx->tp                             = NULL;
    rfc_ctx->tp_cap                         = 0;
    rfc_ctx->tp_cnt                         = 0;
    rfc_ctx->tp_locked                      = 0;
    rfc_ctx->tp_prune_threshold             = (size_t)-1;
    rfc_ctx->tp_prune_size                  = (size_t)-1;
#endif /*RFC_TP_SUPPORT*/


#if RFC_HCM_SUPPORT
    if( rfc_ctx->class_count )
    {
        /* HCM method initialization */
        rfc_ctx->internal.hcm.IZ            = 0;
        rfc_ctx->internal.hcm.IR            = 1;
        /* Residue */
        rfc_ctx->internal.hcm.stack_cap     = 2 * rfc_ctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */
        rfc_ctx->internal.hcm.stack         = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, rfc_ctx->internal.hcm.stack_cap, 
                                                                                      sizeof(rfc_value_tuple_s), RFC_MEM_AIM_HCM );
    }
#endif /*RFC_HCM_SUPPORT*/

#if RFC_AT_SUPPORT
    rfc_ctx->at.Sa                          = NULL;
    rfc_ctx->at.Sm                          = NULL;
    rfc_ctx->at.count                       = 0;
    rfc_ctx->at.M                           = 0.0;
    rfc_ctx->at.Sm_rig                      = 0.0;
    rfc_ctx->at.R_rig                       = 0.0;
    rfc_ctx->at.R_pinned                    = false;

    rfc_ctx->internal.at_haigh.count        = 0;
#endif /*RFC_AT_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT;

#if RFC_DAMAGE_FAST
    if( rfc_ctx->class_count )
    {
        rfc_ctx->damage_lut                 = (double*)rfc_ctx->mem_alloc( rfc_ctx->damage_lut,    class_count * class_count, 
                                                                           sizeof(double), RFC_MEM_AIM_DLUT );
        rfc_ctx->damage_lut_inapt           = 1;
#if RFC_AT_SUPPORT
        rfc_ctx->amplitude_lut              = (double*)rfc_ctx->mem_alloc( rfc_ctx->amplitude_lut, class_count * class_count, 
                                                                           sizeof(double), RFC_MEM_AIM_DLUT );
#endif /*RFC_AT_SUPPORT*/
        damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}


/**
 * @brief      Initialize Woehler parameters to Miners' elementary rule
 *
 * @param      ctx   The rfc context
 * @param      sx    The amplitude "SA"
 * @param      nx    The cycles "N" according to Sa
 * @param      k     The slope "k"
 *
 * @return     true on success
 */
bool RFC_wl_init_elementary( void *ctx, double sx, double nx, double k )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        rfc_ctx->damage_lut_inapt++;
    }
#endif /*RFC_DAMAGE_FAST*/

/* Woehler curve */
    rfc_ctx->wl_sx        =  sx;
    rfc_ctx->wl_nx        =  nx;
    rfc_ctx->wl_k         = -fabs(k);
#if !RFC_MINIMAL
    rfc_ctx->wl_sd        =  0.0;            /* No fatigue strength */
    rfc_ctx->wl_nd        =  DBL_MAX;
    rfc_ctx->wl_k2        =  rfc_ctx->wl_k;
    rfc_ctx->wl_q         =  fabs(k) - 1;    /* Default value for fatigue strength depression */
    rfc_ctx->wl_q2        =  rfc_ctx->wl_q;  /* Default value for fatigue strength depression */
    rfc_ctx->wl_omission  =  0.0;            /* No omission per default */

    /* Make a shadow copy of the Woehler parameters */
    RFC_wl_param_get( rfc_ctx, &rfc_ctx->internal.wl );
#endif /*!RFC_MINIMAL*/

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}


#if !RFC_MINIMAL
/**
 * @brief      Initialize Woehler parameters to Miners' original rule
 *
 * @param      ctx   The rfc context
 * @param      sd    The amplitude "SD"
 * @param      nd    The cycles "ND"
 * @param      k     The slope "k"
 *
 * @return     true on success
 */
bool RFC_wl_init_original( void *ctx, double sd, double nd, double k )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !RFC_wl_init_elementary( ctx, sd, nd, k ) )
    {
        return false;
    }

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        rfc_ctx->damage_lut_inapt++;
    }
#endif /*RFC_DAMAGE_FAST*/

/* Woehler curve */
    rfc_ctx->wl_sd = sd;             /* Points sx/nx and sd/nd coincide */
    rfc_ctx->wl_nd = nd;

    /* Make a shadow copy of the Woehler parameters */
    RFC_wl_param_get( rfc_ctx, &rfc_ctx->internal.wl );

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}


/**
 * @brief      Initialize Woehler parameters to Miners' modified rule
 *
 * @param      ctx   The rfc context
 * @param      sx    The amplitude "Sa"
 * @param      nx    The cycles "N" according to Sa
 * @param      k     The slope before sx/nx
 * @param      k2    The slope after sx/nx
 *
 * @return     true on success
 */
bool RFC_wl_init_modified( void *ctx, double sx, double nx, double k, double k2 )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !RFC_wl_init_elementary( ctx, sx, nx, k ) )
    {
        return false;
    }

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        rfc_ctx->damage_lut_inapt++;
    }
#endif /*RFC_DAMAGE_FAST*/

/* Woehler curve */
    rfc_ctx->wl_k2 =  rfc_ctx->wl_k2;   /* Woehler slope after sx/nx */
    rfc_ctx->wl_q2 =  fabs(k2) - 1;     /* Default value for fatigue strength depression */

    /* Make a shadow copy of the Woehler parameters */
    RFC_wl_param_get( rfc_ctx, &rfc_ctx->internal.wl );

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}


/**
 * @brief      Initialize Woehler curve
 *
 * @param      ctx       The rainflow context
 * @param[in]  wl_param  The parameter set
 *
 * @return     true on success
 */
bool RFC_wl_init_any( void *ctx, const rfc_wl_param_s* wl_param )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !wl_param || !RFC_wl_init_elementary( ctx, wl_param->sx, wl_param->nx, wl_param->k ) )
    {
        return false;
    }

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        rfc_ctx->damage_lut_inapt++;
    }
#endif /*RFC_DAMAGE_FAST*/

/* Woehler curve */
    rfc_ctx->wl_sd       =  wl_param->sd;
    rfc_ctx->wl_nd       =  wl_param->nd;
    rfc_ctx->wl_k2       = -fabs( wl_param->k2 );      /* "Miner elementary", if k == k2 */
    rfc_ctx->wl_q2       =  fabs( wl_param->k2 ) - 1;
    rfc_ctx->wl_omission =  wl_param->omission;

    /* Make a shadow copy of the Woehler parameters */
    RFC_wl_param_get( rfc_ctx, &rfc_ctx->internal.wl );

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}
#endif /*!RFC_MINIMAL*/


#if RFC_TP_SUPPORT
/**
 * @brief      Initialize tp buffer
 *
 * @param      ctx        The rainflow context
 * @param      tp         The buffer for tp
 * @param      tp_cap     The tp capability
 * @param      is_static  Indicates if tp is static
 *
 * @return     true on success
 */
bool RFC_tp_init( void *ctx, rfc_value_tuple_s *tp, size_t tp_cap, bool is_static )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( !tp )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->tp     = tp;
    rfc_ctx->tp_cap = tp_cap;
    rfc_ctx->tp_cnt = 0;
    
    rfc_ctx->internal.tp_static = is_static;

    return true;
}


/**
 * @brief      Initialize autoprune parameters
 *
 * @param      ctx        The rainflow context
 * @param      autoprune  The flag for autopruning
 * @param      size       The size to prune to
 * @param      threshold  The threshold when to prune
 *
 * @return     true on success
 */
bool RFC_tp_init_autoprune( void *ctx, bool autoprune, size_t size, size_t threshold )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->internal.flags      = rfc_ctx->internal.flags & ~RFC_FLAGS_TPAUTOPRUNE | (autoprune ? RFC_FLAGS_TPAUTOPRUNE : 0);
    rfc_ctx->tp_prune_threshold  = threshold;
    rfc_ctx->tp_prune_size       = size;

    return true;
}


/**
 * @brief      Drop turning points from storage, to avoid memory excess
 *
 * @param      ctx    The rainflow context
 * @param      limit  The excepted number of points left in turning points
 *                    storage (May be more, if residuals aren't neglected)
 * @param      flags  The flags (see RFC_FLAGS_TPPRUNE_...)
 *
 * @return     true on success
 */
bool RFC_tp_prune( void *ctx, size_t limit, rfc_flags_e flags )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->tp_cnt > limit )
    {
        rfc_value_tuple_s   *src_beg_it,    /* Source (begin) (tp) */
                            *src_end_it,    /* Source (end) (tp) */
                            *src_it,        /* Source iterator (tp) */
                            *dst_it,        /* Destination iterator (tp) */
                            *res_it;        /* Residue iterator (res) */
        size_t               src_i,         /* Source, position in tp base 0 */
                             dst_i,         /* New turning points, index base 0 (tp) */
                             res_i;         /* Residue, index base 0 (res) */

        size_t               removal;       /* Number of turning points to remove */
        size_t               pos_offset;    /* First position (stream) in tp and dh, base 0 */
        bool                 preserve_pos;  /* Don't justify position */
        bool                 preserve_res;  /* Don't remove turning points, if referenced by residue */

        removal     = rfc_ctx->tp_cnt - limit;
        dst_it      = rfc_ctx->tp;
        dst_i       = 0;
        src_beg_it  = rfc_ctx->tp + removal;
        src_end_it  = rfc_ctx->tp + rfc_ctx->tp_cnt
                      + ( ( rfc_ctx->state == RFC_STATE_BUSY_INTERIM ) ? 1 : 0 );
        src_it      = src_beg_it;
        src_i       = removal;
        res_it      = rfc_ctx->residue;
        res_i       = 0;
        pos_offset  = 0;

        preserve_pos = ( flags & RFC_FLAGS_TPPRUNE_PRESERVE_POS ) > 0;  /* Preserve (stream) position */
        preserve_res = ( flags & RFC_FLAGS_TPPRUNE_PRESERVE_RES ) > 0;  /* Preserve residual turning points */

        /* Move turning points ahead */
        while( src_it < src_end_it || res_i < rfc_ctx->residue_cnt )
        {
            /* Check if there are still residual points to consider */
            while( res_i < rfc_ctx->residue_cnt && res_it->tp_pos <= src_i + 1 )
            {
                /* Check if residue refers a turning point from removal area */

                /* First new turning point delivers new offset */
                if( !res_i && !preserve_pos )
                {
                    pos_offset = res_it->pos;  // pos is base 1
                    assert( pos_offset );
                    pos_offset--;
                }

                /* Residual point refers current source position? */
                if( res_it->tp_pos == src_i + 1 )
                {
                    /* Turning will be processed in this inner loop */
                    src_it++;
                    src_i++;
                }

                if( preserve_res )
                {
                    /* Adjust residue reference information */
                    res_it->pos -= pos_offset;

                    /* Set residual turning point */
                    if( !tp_set( rfc_ctx, dst_i + 1, res_it ) )
                    {
                        return error_raise( rfc_ctx, RFC_ERROR_TP );
                    }

                    dst_it++;
                    dst_i++;
                    res_it++;
                    res_i++;
                }
                else
                {
                    /* Residual turning point refers first point now */
                    res_it->tp_pos  = 0;  /* Index 0 => "none" */
                    res_it->pos    -= pos_offset;
                    res_it++;
                    res_i++;
                }
            }

            if( src_it < src_end_it )
            {
                rfc_value_tuple_s *cpy;

                /* First new turning point delivers new offset */
                if( !dst_i && !preserve_pos )
                {
                    pos_offset = src_it->pos;
                    assert( pos_offset );
                    pos_offset--;
                }

                /* Copy turning point from source */
                if( !tp_get( rfc_ctx, src_i + 1, &cpy ) )
                {
                    return error_raise( rfc_ctx, RFC_ERROR_TP );
                }

                /* Adjust stream position */
                cpy->pos -= pos_offset;

                /* Move turning point in tp stack */
                if( !tp_set( rfc_ctx, dst_i + 1, cpy ) )
                {
                    return error_raise( rfc_ctx, RFC_ERROR_TP );
                }

                dst_it++;
                dst_i++;
                src_it++;
                src_i++;
            }
        }

        rfc_ctx->tp_cnt                  = dst_i;
        rfc_ctx->internal.pos           -= pos_offset;
        rfc_ctx->internal.pos_offset    += pos_offset;

#if RFC_DH_SUPPORT
        /* Shift damage history */
        if( rfc_ctx->dh && pos_offset )
        {
            assert( rfc_ctx->dh_cnt >= pos_offset );
            memcpy( rfc_ctx->dh, rfc_ctx->dh + pos_offset, rfc_ctx->dh_cnt - pos_offset );
            rfc_ctx->dh_cnt -= pos_offset;
        }
#endif /*RFC_DH_SUPPORT*/
    }
    
    return true;
}


/**
 * @brief      Clear turning point storage
 *
 * @param      ctx   The rainflow context
 *
 * @return     true on success
 */
bool RFC_tp_clear( void *ctx )
{
    size_t i;

    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->tp_cnt = 0;

    for( i = 0; i < rfc_ctx->residue_cnt; i++ )
    {
        rfc_ctx->residue[i].tp_pos = 0;
    }

    return true;
}

#endif /*RFC_TP_SUPPORT*/


#if RFC_DH_SUPPORT
/**
 * @brief      Initialize damage history storage
 *
 * @param      ctx        The rainflow context
 * @param[in]  method     The mode, how to spread (RFC_SD_...)
 * @param      dh         The storage buffer
 * @param      dh_cap     The capacity of dh
 * @param      is_static  true, if dh is static and should not be freed
 *
 * @return     true on success
 */
bool RFC_dh_init( void *ctx, rfc_sd_method_e method, double *dh, size_t dh_cap, bool is_static )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->spread_damage_method = method;
    rfc_ctx->dh                   = dh;
    rfc_ctx->dh_cap               = dh_cap;
    rfc_ctx->dh_cnt               = 0;

    rfc_ctx->internal.dh_static   = is_static;

    return true;
}
#endif /*RFC_DH_SUPPORT*/


#if RFC_AT_SUPPORT
/**
 * @brief      Initialize amplitude transformation
 *
 * @param      ctx        The rainflow context
 * @param      Sa         The reference curve vector, amplitude part
 * @param      Sm         The reference curve vector, mean load part. If Sa and
 *                        Sm_norm are NULL, the standard (FKM) is applied
 * @param      count      The capacity of Sa and Sm
 * @param      M          The mean stress sensitivity
 * @param      Sm_rig     The mean load applied on the test rig
 * @param      R_rig      The mean load ratio applied on the test rig
 * @param      R_pinned   true, if R is constant on test rig (R_rig is used).
 *                        false if Sm is constant on test rig (Sm_rig is used)
 * @param      symmetric  true if Haigh diagram is symmetric at Sa(R=-1)
 *
 * @return     true on success
 */
bool RFC_at_init( void *ctx, const double *Sa, const double *Sm, unsigned count, 
                             double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( M < 0.0 )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    if( count )
    {
        /* Reference curve given, doing some checks */

        unsigned n;

        if( !Sa || !Sm || symmetric || count < 2 )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }

        /* Check for valid input */
        for( n = 0; n < count; n++ )
        {
            if( Sa[n] <= 0.0 ) break;

            if( !n ) continue;

            if( Sm[n-1] >= Sm[n] || Sm[n-1] / Sa[n-1] > Sm[n] / Sa[n] ) break;
        }

        if( n < count )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }

        rfc_ctx->at.Sa       = Sa;
        rfc_ctx->at.Sm       = Sm;
        rfc_ctx->at.count    = count;
        rfc_ctx->at.M        = M;
        rfc_ctx->at.Sm_rig   = Sm_rig;
        rfc_ctx->at.R_rig    = R_rig;
        rfc_ctx->at.R_pinned = R_pinned;
    }
    else
    {
        /* No reference curve given */

        double Sa_R_Inf, Sa_R_0, Sa_R_0p5;

        assert( !Sa && !Sm );

        if( M > 0.0 )
        {
            Sa_R_Inf = 1.0 / ( 1.0 - M );                      /* y = -x && y = Sa(R=-1) - Mx                  */
            Sa_R_0   = 1.0 / ( 1.0 + M );                      /* y =  x && y = Sa(R=-1) - Mx                  */
            Sa_R_0p5 = Sa_R_0 * ( 1.0 + M/3 ) / ( 1.0 + M );   /* Backtrace Sa(R=0) to Sa(R=-1) with M/3, then */
                                                               /* 3y = x && y = Sa(R=-1) - (M/3)x              */
            if( symmetric )
            {
                /* Build symmetrical reference curve */
                /* Symmetric around R=-1 (Sm=0) */

                double *Sa_ = rfc_ctx->internal.at_haigh.Sa;
                double *Sm_ = rfc_ctx->internal.at_haigh.Sm;

                assert( NUMEL( rfc_ctx->internal.at_haigh.Sa ) >= 5 );
                
                rfc_ctx->internal.at_haigh.count = 5;

                Sa_[0] = Sa_R_0p5; Sm_[0] = -Sa_R_0p5 * 3.0;
                Sa_[1] = Sa_R_0;   Sm_[1] = -Sa_R_0;
                Sa_[2] = 1.0;      Sm_[2] =  0.0;
                Sa_[3] = Sa_[1];   Sm_[3] = -Sm_[1];
                Sa_[4] = Sa_[0];   Sm_[4] = -Sm_[0];
            }
            else
            {
                /* Build non-symmetric reference curve */
                double *Sa_ = rfc_ctx->internal.at_haigh.Sa;
                double *Sm_ = rfc_ctx->internal.at_haigh.Sm;

                assert( NUMEL( rfc_ctx->internal.at_haigh.Sa ) >= 3 );
                
                rfc_ctx->internal.at_haigh.count = 3;

                Sa_[0] = Sa_R_Inf; Sm_[0] = -Sa_R_Inf;
                Sa_[1] = Sa_R_0;   Sm_[1] =  Sa_R_0;
                Sa_[2] = Sa_R_0p5; Sm_[2] =  Sa_R_0p5 * 3.0;
            }

            rfc_ctx->at.Sa       = rfc_ctx->internal.at_haigh.Sa;
            rfc_ctx->at.Sm       = rfc_ctx->internal.at_haigh.Sm;
            rfc_ctx->at.count    = rfc_ctx->internal.at_haigh.count;
        }
        else
        {
            rfc_ctx->at.Sa       = NULL;
            rfc_ctx->at.Sm       = NULL;
            rfc_ctx->at.count    = 0;
        }

        rfc_ctx->at.M        = M;
        rfc_ctx->at.Sm_rig   = Sm_rig;
        rfc_ctx->at.R_rig    = R_rig;
        rfc_ctx->at.R_pinned = R_pinned;
    }

#if RFC_DAMAGE_FAST
    damage_lut_init( rfc_ctx );
#endif /*RFC_DAMAGE_FAST*/

    return true;
}
#endif /*RFC_AT_SUPPORT*/


#if !RFC_MINIMAL
/**
 * @brief      Clear all data generated while counting
 *
 * @param      rfc_ctx  The rainflow context
 */
bool RFC_clear_counts( void *ctx )
{
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    if( rfc_ctx->rfm )
    {
        memset( rfc_ctx->rfm, 0, sizeof(rfc_counts_t) * rfc_ctx->class_count * rfc_ctx->class_count );
    }

    if( rfc_ctx->rp )
    {
        memset( rfc_ctx->rp, 0, sizeof(rfc_counts_t) * rfc_ctx->class_count );
    }

    if( rfc_ctx->lc )
    {
        memset( rfc_ctx->lc, 0, sizeof(rfc_counts_t) * rfc_ctx->class_count );
    }

    rfc_ctx->residue_cnt                = 0;

    rfc_ctx->internal.slope             = 0;
    rfc_ctx->internal.extrema[0]        = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]        = nil;  /* local maximum */
#if RFC_GLOBAL_EXTREMA
    rfc_ctx->internal.extrema_changed   = false;
#endif
    rfc_ctx->internal.pos               = 0;
    rfc_ctx->internal.pos_offset        = 0;
    
    rfc_ctx->damage                     = 0.0;
    rfc_ctx->damage_residue             = 0.0;

#if RFC_HCM_SUPPORT
    /* Reset stack pointers */
    rfc_ctx->internal.hcm.IR            = 1;
    rfc_ctx->internal.hcm.IZ            = 0;
#endif /*RFC_HCM_SUPPORT*/

#if RFC_TP_SUPPORT
    /* rfc_ctx->tp_cnt is set to zero, but turning points are still available */
    rfc_ctx->internal.margin[0]         = nil;  /* left margin */
    rfc_ctx->internal.margin[1]         = nil;  /* right margin */
    rfc_ctx->internal.margin_stage      = 0;
    rfc_ctx->tp_cnt                     = 0;
    rfc_ctx->tp_locked                  = 0;
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    rfc_ctx->dh_cnt                     = 0;
#endif /*RFC_DH_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT;

    return true;
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      De-initialization (freeing memory).
 *
 * @param      ctx   The rainflow context
 *
 * @return     true on success
 */
bool RFC_deinit( void *ctx )
{
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    if( !rfc_ctx->internal.res_static &&
        rfc_ctx->residue )              rfc_ctx->mem_alloc( rfc_ctx->residue,       0, 0, RFC_MEM_AIM_RESIDUE );
    if( rfc_ctx->rfm )                  rfc_ctx->mem_alloc( rfc_ctx->rfm,           0, 0, RFC_MEM_AIM_MATRIX );
#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )           rfc_ctx->mem_alloc( rfc_ctx->damage_lut,    0, 0, RFC_MEM_AIM_DLUT );
#if RFC_AT_SUPPORT
    if( rfc_ctx->amplitude_lut )        rfc_ctx->mem_alloc( rfc_ctx->amplitude_lut, 0, 0, RFC_MEM_AIM_DLUT );
#endif /*RFC_AT_SUPPORT*/
#endif /*RFC_DAMAGE_FAST*/
#if !RFC_MINIMAL
    if( rfc_ctx->rp )                   rfc_ctx->mem_alloc( rfc_ctx->rp,            0, 0, RFC_MEM_AIM_RP );
    if( rfc_ctx->lc )                   rfc_ctx->mem_alloc( rfc_ctx->lc,            0, 0, RFC_MEM_AIM_LC );
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
    if( rfc_ctx->tp && !rfc_ctx->internal.tp_static )
    {
                                        rfc_ctx->mem_alloc( rfc_ctx->tp,            0, 0, RFC_MEM_AIM_TP );
    }           
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
    if( rfc_ctx->dh && !rfc_ctx->internal.dh_static )
    {               
                                        rfc_ctx->mem_alloc( rfc_ctx->dh,            0, 0, RFC_MEM_AIM_DH );
    }
#endif /*RFC_DH_SUPPORT*/

#if RFC_DAMAGE_FAST
    rfc_ctx->damage_lut                 = NULL;
    rfc_ctx->damage_lut_inapt           = 1;
#if RFC_AT_SUPPORT
    rfc_ctx->amplitude_lut              = NULL;
#endif /*RFC_AT_SUPPORT*/
#endif /*RFC_DAMAGE_FAST*/

    rfc_ctx->residue                    = NULL;
    rfc_ctx->residue_cap                = 0;
    rfc_ctx->residue_cnt                = 0;

    rfc_ctx->rfm                        = NULL;
#if !RFC_MINIMAL
    rfc_ctx->rp                         = NULL;
    rfc_ctx->lc                         = NULL;
#endif /*!RFC_MINIMAL*/
    
    rfc_ctx->internal.slope             = 0;
    rfc_ctx->internal.extrema[0]        = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]        = nil;  /* local maximum */
#if RFC_GLOBAL_EXTREMA
    rfc_ctx->internal.extrema_changed   = false;
#endif /*RFC_GLOBAL_EXTREMA*/
    rfc_ctx->internal.pos               = 0;
    rfc_ctx->internal.pos_offset        = 0;
#if RFC_TP_SUPPORT
    rfc_ctx->internal.margin[0]         = nil;  /* left margin */
    rfc_ctx->internal.margin[1]         = nil;  /* right margin */
    rfc_ctx->internal.margin_stage      = 0;

    rfc_ctx->tp                         = NULL;
    rfc_ctx->tp_cap                     = 0;
    rfc_ctx->tp_cnt                     = 0;
    rfc_ctx->tp_locked                  = 0;
    rfc_ctx->internal.tp_static         = false;
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    rfc_ctx->dh                         = NULL;
    rfc_ctx->dh_cap                     = 0;
    rfc_ctx->dh_cnt                     = 0;
    rfc_ctx->internal.dh_static         = false;
#endif /*RFC_DH_SUPPORT*/

#if RFC_AT_SUPPORT
    rfc_ctx->at.Sa                      = NULL;
    rfc_ctx->at.Sm                      = NULL;
    rfc_ctx->at.count                   = 0;
    rfc_ctx->at.M                       = 0.0;
    rfc_ctx->at.Sm_rig                  = 0.0;
    rfc_ctx->at.R_rig                   = 0.0;
    rfc_ctx->at.R_pinned                = false;

    rfc_ctx->internal.at_haigh.count    = 0;
#endif /*RFC_AT_SUPPORT*/

#if RFC_HCM_SUPPORT
    /* Remove stack */
    if( rfc_ctx->internal.hcm.stack )   rfc_ctx->mem_alloc( rfc_ctx->internal.hcm.stack, 0, 0, RFC_MEM_AIM_HCM );

    rfc_ctx->internal.hcm.stack         = NULL;
    rfc_ctx->internal.hcm.stack_cap     = 0;

    /* Stack pointers */
    rfc_ctx->internal.hcm.IZ            = 0;
    rfc_ctx->internal.hcm.IR            = 1;

#endif /*RFC_HCM_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT0;

    return true;
}


/**
 * @brief      "Feed" counting algorithm with data samples (consecutive calls
 *             allowed).
 *
 * @param      ctx         The rainflow context
 * @param[in]  data        The data
 * @param      data_count  The data count
 *
 * @return     true on success
 */
bool RFC_feed( void *ctx, const rfc_value_t * data, size_t data_count )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( data_count && !data ) return false;

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    /* Process data */
    while( data_count-- )
    {
        rfc_value_tuple_s tp = { *data++ };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

        /* Assign class and global position (base 1) */
        tp.pos = ++rfc_ctx->internal.pos;
        tp.cls = QUANTIZE( rfc_ctx, tp.value );

        if( tp.cls >= rfc_ctx->class_count && rfc_ctx->class_count )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
        
        if( !feed_once( rfc_ctx, &tp, rfc_ctx->internal.flags ) ) return false;
    }

    return true;
}


/**
 * @brief      Do countings for a given cycle
 *
 * @param      ctx       The rainflow context
 * @param[in]  from_val  The from value
 * @param[in]  to_val    The to value
 * @param[in]  flags     The flags
 *
 * @return     true on success
 */
#if !RFC_MINIMAL
bool RFC_cycle_process_counts( void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_flags_e flags )
{
    rfc_value_tuple_s from = {from_val}, to = {to_val};
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    from.cls = QUANTIZE( rfc_ctx, from_val );
    to.cls   = QUANTIZE( rfc_ctx, to_val );

    cycle_process_counts( rfc_ctx, &from, &to, /*next*/ NULL, flags );

    return true;
}


/**
 * @brief      "Feed" counting algorithm with data samples, scaled by a factor
 *             (consecutive calls allowed).
 *
 * @param      ctx         The rainflow context
 * @param[in]  data        The data
 * @param      data_count  The data count
 * @param      factor      The factor
 *
 * @return     true on success
 */
bool RFC_feed_scaled( void *ctx, const rfc_value_t * data, size_t data_count, double factor )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( data_count && !data ) return false;

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    /* Process data */
    while( data_count-- )
    {
        rfc_value_tuple_s tp = { *data++ * factor };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

        /* Assign class and global position (base 1) */
        tp.cls = QUANTIZE( rfc_ctx, tp.value );
        tp.pos = ++rfc_ctx->internal.pos;

        if( tp.cls >= rfc_ctx->class_count && rfc_ctx->class_count )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
        
        if( !feed_once( rfc_ctx, &tp, rfc_ctx->internal.flags ) ) return false;
    }

    return true;
}


/**
 * @brief         Feed counting algorithm with data tuples.
 *
 * @param         ctx         The rainflow context
 * @param[in,out] data        The data tuples
 * @param         data_count  The data count
 *
 * @return        true on success
 */
bool RFC_feed_tuple( void *ctx, rfc_value_tuple_s *data, size_t data_count )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( data_count && !data ) return false;

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    /* Process data */
    while( data_count-- )
    {
        if( data->cls >= rfc_ctx->class_count && rfc_ctx->class_count )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
        
        if( !feed_once( rfc_ctx, data++, rfc_ctx->internal.flags ) ) return false;
    }

    return true;
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Finalize pending counts and turning point storage.
 *
 * @param      ctx              The rainflow context
 * @param      residual_method  The residual method (RFC_RES_...)
 *
 * @return     true on success
 */
bool RFC_finalize( void *ctx, rfc_res_method_e residual_method )
{
    double damage;
    bool ok;
    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    damage = rfc_ctx->damage;

#if RFC_USE_DELEGATES
    if( rfc_ctx->finalize_fcn )
    {
        ok = rfc_ctx->finalize_fcn( rfc_ctx, residual_method );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        int flags = rfc_ctx->internal.flags;

#if !RFC_MINIMAL
        /* Level crossing counting is already considered for residue */
        flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

        switch( residual_method )
        {
            case RFC_RES_NONE:
                /* FALLTHROUGH */
            case RFC_RES_IGNORE:
                ok = finalize_res_ignore( rfc_ctx, flags );
                break;
#if !RFC_MINIMAL
            case RFC_RES_DISCARD:
                ok = finalize_res_discard( rfc_ctx, flags );
                break;
            case RFC_RES_HALFCYCLES:
                ok = finalize_res_weight_cycles( rfc_ctx, rfc_ctx->half_inc, flags );
                break;
            case RFC_RES_FULLCYCLES:
                ok = finalize_res_weight_cycles( rfc_ctx, rfc_ctx->full_inc, flags );
                break;
            case RFC_RES_CLORMANN_SEEGER:
                ok = finalize_res_clormann_seeger( rfc_ctx, flags );
                break;
            case RFC_RES_REPEATED:
                ok = finalize_res_repeated( rfc_ctx, flags );
                break;
            case RFC_RES_RP_DIN45667:
                ok = RFC_finalize_res_rp_DIN45667( rfc_ctx, flags );
                break;
#endif /*!RFC_MINIMAL*/
            default:
                assert( false );
                ok = error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
        assert( rfc_ctx->state == RFC_STATE_FINALIZE );
    }

#if !RFC_MINIMAL
    if( rfc_ctx->counting_method == RFC_COUNTING_METHOD_NONE || !rfc_ctx->class_count )
    {
#else /*RFC_MINIMAL*/
    if( !rfc_ctx->class_count )
    {
#endif /*!RFC_MINIMAL*/
        rfc_ctx->residue_cnt = 0;
    }

    rfc_ctx->damage_residue = rfc_ctx->damage - damage;
    rfc_ctx->state          = ok ? RFC_STATE_FINISHED : RFC_STATE_ERROR;

    return ok;
}


#if !RFC_MINIMAL
/**
 * @brief      Make rainflow matrix symmetrical
 *
 * @param      ctx   The rainflow context
 *
 * @return     true on success
 */
bool RFC_rfm_make_symmetric( void *ctx )
{
    unsigned       class_count;
    unsigned       from, to;
    rfc_counts_t  *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    rfm = rfc_ctx->rfm;

    if( !rfm )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    /* Cumulate entries symmetric by major diagonal */
    for( from = 0; from < class_count; from++ )
    {
        for( to = from + 1; to < class_count; to++ )
        {
            // to > from, always
            rfm[ class_count * from + to   ] += rfm[ class_count * to + from ];
            rfm[ class_count * to   + from ]  = 0;
        }
    }

    return true;
}


/**
 * @brief      Returns the number of non zero entries in rainflow matrix
 *
 * @param[in]  ctx   The rainflow context
 *
 * @return     true on success
 */
bool RFC_rfm_non_zeros( const void *ctx, unsigned *count )
{
    unsigned            class_count;
    unsigned            from, to;
    rfc_counts_t       *rfm_it;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !count )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm_it = rfc_ctx->rfm;

    if( !rfm_it || !class_count )
    {
        return false;
    }

    *count     = 0;
    for( from = 0; from < class_count; from++ )
    {
        for( to = 0; to < class_count; to++, rfm_it )
        {
            if( *rfm_it ) *count++;
        }
    }

    return true;
}


/**
 * @brief      Get the rainflow matrix as sparse elements
 *
 * @param      ctx     The rainflow context
 * @param[out] buffer  The elements buffer, if NULL memory will be allocated
 * @param[out] count   The number of elements in buffer
 *
 * @return     true on success
 * @note       The counts are natively returned, regardless of .full_inc!
*/
bool RFC_rfm_get( const void *ctx, rfc_rfm_item_s **buffer, unsigned *count )
{
    unsigned            class_count;
    unsigned            from, to;
    unsigned            count_old;
    rfc_counts_t       *rfm_it;
    rfc_rfm_item_s     *item;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !buffer || !count )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm_it = rfc_ctx->rfm;

    if( !rfm_it || !class_count )
    {
        return false;
    }

    // *buffer = NULL;
    count_old  = *count;
    *count     = 0;
    for( from = 0; from < class_count; from++ )
    {
        for( to = 0; to < class_count; to++, rfm_it++ )
        {
            if( *rfm_it )
            {
                (*count)++;
            }
        }
    }

    if( *count > count_old )
    {
        *buffer = rfc_ctx->mem_alloc( *buffer, *count, sizeof(rfc_rfm_item_s), RFC_MEM_AIM_RFM_ELEMENTS );

        if( !*buffer )
        {
            error_raise( rfc_ctx, RFC_ERROR_MEMORY );
            return false;
        }
    }
        
    item   = *buffer;
    rfm_it = rfc_ctx->rfm;
    for( from = 0; from < class_count; from++ )
    {
        for( to = 0; to < class_count; to++, rfm_it++ )
        {
            if( *rfm_it )
            {
                item->from   = from;
                item->to     = to;
                item->counts = *rfm_it;

                item++;
            }
        }
    }

    return true;
}


/**
 * @brief      Set (or increment) rainflow matrix with given elements
 *
 * @param      ctx       The rainflow context
 * @param[in]  buffer    The elements buffer
 * @param      count     The number of elements in buffer
 * @param      add_only  Counts are added if set to true
 *
 * @return     true on success
 * @note       The counts are natively added, regardless of .full_inc!
 */
bool RFC_rfm_set( void *ctx, const rfc_rfm_item_s *buffer, unsigned count, bool add_only )
{
          unsigned           class_count, i;
    const rfc_rfm_item_s    *item;
          rfc_counts_t      *rfm;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !buffer )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    if( !add_only )
    {
        memset( rfm, (rfc_counts_t)0, sizeof(rfc_rfm_item_s) * class_count * class_count );
    }

    item = buffer;
    for( i = 0; i < count; i++ )
    {
        unsigned from, to;

        assert( item->from >= rfc_ctx->class_offset );
        assert( item->to   >= rfc_ctx->class_offset );

        from = QUANTIZE( rfc_ctx, item->from );
        to   = QUANTIZE( rfc_ctx, item->to );

        if( from > class_count ) from = class_count;
        if( to   > class_count ) to   = class_count;

        rfm[ from * class_count + to ] += item->counts;
    }

    return true;
}


/**
 * @brief      Get counts of a single element from the rainflow matrix
 *
 * @param      ctx       The rainflow context
 * @param      from_val  The cycles start value
 * @param      to_val    The cycles target value
 * @param[out] counts    The corresponding count from the matrix element (not cycles!), regardless of .full_inc!
 *
 * @return     true on success
 */
bool RFC_rfm_peek( const void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t *counts )
{
    unsigned           from, to;
    unsigned           class_count;
    rfc_counts_t      *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    assert( from_val >= rfc_ctx->class_offset );
    assert( to_val   >= rfc_ctx->class_offset );

    from = QUANTIZE( rfc_ctx, from_val );
    to   = QUANTIZE( rfc_ctx, to_val );

    if( from > class_count ) from = class_count;
    if( to   > class_count ) to   = class_count;

    if( counts )
    {
        *counts = rfm[ from * class_count + to ];
    }

    return true;
}


/**
 * @brief      Set (or increment) one matrix value of the rainflow matrix
 *
 * @param      ctx       The rainflow context
 * @param      from_val  The cycles start value
 * @param      to_val    The cycles target value
 * @param      counts    The count value for the matrix element (not cycles!), regardless of .full_inc!
 * @param      add_only  Value is added if set to true
 *
 * @return     true on success
 */
bool RFC_rfm_poke( void *ctx, rfc_value_t from_val, rfc_value_t to_val, rfc_counts_t counts, bool add_only )
{
    unsigned           from, to;
    unsigned           class_count;
    rfc_counts_t      *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    assert( from_val >= rfc_ctx->class_offset );
    assert( to_val   >= rfc_ctx->class_offset );

    from = QUANTIZE( rfc_ctx, from_val );
    to   = QUANTIZE( rfc_ctx, to_val );

    if( from > class_count ) from = class_count;
    if( to   > class_count ) to   = class_count;

    if( add_only )
    {
        rfm[ from * class_count + to ] += counts;
    }
    else
    {
        rfm[ from * class_count + to ] = counts;
    }

    return true;
}


/**
 * @brief      Sum cycles of a rainflow matrix region
 *
 * @param      ctx         The rainflow context
 * @param      from_first  The first start class (row)
 * @param      from_last   The last start class (row)
 * @param      to_first    The first target class (col)
 * @param      to_last     The last target class (col)
 * @param      count       The sum of the matrix region
 *
 * @return     true on success
 * @note       The sum is natively built, regardless of .full_inc!
 */
bool RFC_rfm_sum( const void *ctx, unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, rfc_counts_t *count )
{
    unsigned         from, to;
    unsigned         class_count;
    rfc_counts_t    *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    assert( from_first < class_count );
    assert( from_last  < class_count );
    assert( to_first   < class_count );
    assert( to_last    < class_count );
    assert( from_first < from_last );
    assert( to_first   < to_last );

    if( count )
    {
        rfc_counts_t sum = 0;

        for( from = from_first; from <= from_last; from++ )
        {
            for( to = to_first; to < to_last; to++ )
            {
                sum += rfm[ from * class_count + to ];
            }
        }

        *count = sum;
    }

    return true;
}


/**
 * @brief      Calculates the sum of damages for a rainflow matrix
 *             region
 *
 * @param      ctx         The rainflow context
 * @param      from_first  The first start class (row)
 * @param      from_last   The last start class (row)
 * @param      to_first    The first target class (col)
 * @param      to_last     The last target class (col)
 * @param[out] damage      The result (sum)
 *
 * @return     true on success
 */
bool RFC_rfm_damage( const void *ctx, unsigned from_first, unsigned from_last, unsigned to_first, unsigned to_last, double *damage )
{
    unsigned          from, to;
    unsigned          class_count;
    rfc_counts_t     *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    assert( from_first < class_count );
    assert( from_last  < class_count );
    assert( to_first   < class_count );
    assert( to_last    < class_count );
    assert( from_first < from_last );
    assert( to_first   < to_last );

    if( damage )
    {
        double sum = 0.0;
        for( from = from_first; from <= from_last; from++ )
        {
            for( to = to_first; to < to_last; to++ )
            {
                rfc_counts_t count = rfm[ from * class_count + to ];
                double damage_i;
                
                if( !damage_calc( rfc_ctx, from, to, &damage_i, NULL /*Sa_ret*/ ) )
                {
                    return false;
                }

                sum += damage_i * count;
            }
        }

        *damage = sum / rfc_ctx->full_inc;
    }

    return true;
}


/**
 * @brief      Check the consistency of the rainflow matrix
 *
 * @param      rfc_ctx  The rfc context
 *
 * @return     true on success
 */
bool RFC_rfm_check( const void *ctx )
{
    unsigned          class_count;
    rfc_counts_t     *rfm;

    RFC_CTX_CHECK_AND_ASSIGN
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }
    else
    {
        int i;

        for( i = 0; i < (int)rfc_ctx->class_count; i++ )
        {
            /* Matrix diagonal must be all zero */
            if( rfc_ctx->rfm[ i * rfc_ctx->class_count + i ] != 0 )
            {
                return false;
            }
        }
    }
    return true;
}


/**
 * @brief      Get level crossing histogram
 *
 * @param      ctx    The rainflow context
 * @param[out] lc     The buffer for LC histogram (counts), .full_inc represents one "count"!
 * @param[out] level  The buffer for LC upper class borders (dropped if NULL)
 *
 * @return     true on success
 */
bool RFC_lc_get( const void *ctx, rfc_counts_t *lc, rfc_value_t *level )
{
    unsigned i;
    unsigned class_count;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !lc )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !rfc_ctx->lc || !lc || !class_count )
    {
        return false;
    }

    for( i = 0; i < class_count; i++ )
    {
        lc[i] = rfc_ctx->lc[i];

        if( level )
        {
            level[i] = CLASS_UPPER( rfc_ctx, i );
        }
    }

    return true;
}


/**
 * @brief      Create level crossing histogram from rainflow matrix
 *
 * @param      ctx     The rainflow context
 * @param[out] lc      The buffer for LC histogram (counts), .full_inc represents one "count"
 * @param[out] level   The buffer for LC upper class borders (dropped if NULL)
 * @param[in]  rfm     The input rainflow matrix, max be NULL
 * @param      flags   The flags
 *
 * @return     true on success
 * @note       Returned lc usually differs from .lc, when counting is finalized with any flag other than RFC_RES_NONE!
 */
bool RFC_lc_from_rfm( const void *ctx, rfc_counts_t *lc, rfc_value_t *level, const rfc_counts_t *rfm, rfc_flags_e flags )
{
    unsigned             from, to, i;
    unsigned             class_count;
    bool                 up = flags & RFC_FLAGS_COUNT_LC_UP;
    bool                 dn = flags & RFC_FLAGS_COUNT_LC_DN;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !lc )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    rfm = rfc_ctx->rfm;

    if( !rfm || !class_count )
    {
        return false;
    }

    for( i = 0; i < class_count; i++ ) 
    {
        rfc_counts_t sum = 0;

        /* First index (0) counts crossings of upper class limit of the first class */
        if( level )
        {
            level[i] = CLASS_UPPER( rfc_ctx, i );
        }
        
        for( from = 0; from <= i; from++ )
        {
            for( to = i + 1; to < class_count; to++ )
            {
                /* One closed cycle has always a rising and a falling slope */
                 
                if( up )
                {
                    /* Count rising slopes */
                    assert( sum < RFC_COUNTS_LIMIT - rfm[ from * class_count + to ] );
                    sum += rfm[ from * class_count + to ];

                    assert( sum < RFC_COUNTS_LIMIT - rfm[ from * class_count + to ] );
                    sum += rfm[ to * class_count + from ];
                }
                if( dn )
                {
                    /* Count falling slopes */
                    assert( sum < RFC_COUNTS_LIMIT - rfm[ from * class_count + to ] );
                    sum += rfm[ from * class_count + to ];

                    assert( sum < RFC_COUNTS_LIMIT - rfm[ from * class_count + to ] );
                    sum += rfm[ to * class_count + from ];
                }
            }
        }

        lc[i] = sum;
    }

    return true;
}


/**
 * Calculate level crossing counts from rainflow matrix, write results to lc
 * histogram buffer.
 *
 * @param      ctx    The rainflow context
 * @param[out] lc     The buffer for LC histogram (counts), .full_inc represents one "count"!
 * @param[out] level  The buffer for LC upper class borders (dropped if NULL)
 * @param      flags  The flags
 *
 * @return     true on success
 */
bool RFC_lc_from_residue( const void *ctx, rfc_counts_t* lc, rfc_value_t *level, rfc_flags_e flags )
{
          unsigned           i;
          unsigned           class_count;
          bool               up   = flags & RFC_FLAGS_COUNT_LC_UP;
          bool               dn   = flags & RFC_FLAGS_COUNT_LC_DN;
    const rfc_value_tuple_s *from;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !lc )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !class_count )
    {
        return false;
    }

    memset( lc, 0, sizeof(rfc_counts_t) * class_count );

    if( level )
    {
        for( i = 0; i < class_count; i++ )
        {
            level[i] = CLASS_UPPER( rfc_ctx, i );
        }
    }

    from = rfc_ctx->residue;
    for( i = 1; i < rfc_ctx->residue_cnt; i++ ) 
    {
        const rfc_value_tuple_s *to         = from + 1;
              unsigned           class_from = from->cls;
              unsigned           class_to   = to->cls;

        /* Level crossing, count rising and falling slopes.
         * Level crossing histogram (vector storage)
         * Counts class upper bound crossings
         * Class upper bound value = (idx+1) * class_width + class_offset
         */
        if( class_from < class_to && up )
        {
            /* Count rising slopes */
            unsigned idx;
            for( idx = class_from; idx < class_to; idx++ )
            {
                assert( lc[idx] <= RFC_COUNTS_LIMIT );
                lc[idx] += rfc_ctx->full_inc;
            }
        }
        else if( class_to < class_from && dn )
        {
            /* Count falling slopes */
            unsigned idx;
            for( idx = class_to; idx < class_from; idx++ )
            {
                assert( lc[idx] <= RFC_COUNTS_LIMIT );
                lc[idx] += rfc_ctx->full_inc;
            }
        }

        from = to;
    }

    return true;
}


/**
 * @brief      Get range pair histogram
 *
 * @param      ctx          The rainflow context
 * @param[out] rp           The histogram (counts), .full_inc represent one "cycle"!
 * @param[out] Sa           The amplitudes (dropped if NULL)
 *
 * @return     true on success
 */
bool RFC_rp_get( const void *ctx, rfc_counts_t *rp, rfc_value_t *Sa )
{
    unsigned i;
    unsigned class_count;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !rp )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !rfc_ctx->rp || !rp || !class_count )
    {
        return false;
    }

    for( i = 0; i < class_count; i++ )
    {
        rp[i] = rfc_ctx->rp[i];

        if( Sa )
        {
            Sa[i] = rfc_ctx->class_width * i / 2;  /* range / 2 */
        }
    }

    return true;
}


/**
 * @brief      Generate range pair histogram from rainflow matrix
 *
 * @param      ctx          The rainflow context
 * @param[out] rp           The buffer for range pair counts (not cycles!), .full_inc represents one "cycle"!
 * @param[out] Sa           The buffer for amplitudes (dropped if NULL)
 * @param[in]  rfm          The input rainflow matrix (may be NULL)
 *
 * @return     true on success
 */
bool RFC_rp_from_rfm( const void *ctx, rfc_counts_t *rp, rfc_value_t *Sa, const rfc_counts_t *rfm )
{
    unsigned    i, j;
    unsigned    class_count;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !rp )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !rfm )
    {
        rfm = rfc_ctx->rfm;
    }

    if( !rfm || !class_count )
    {
        return false;
    }

    for( i = 0; i < class_count; i++ ) 
    {
        rfc_counts_t sum = (rfc_counts_t)0;

        if( Sa )
        {
            Sa[i] = rfc_ctx->class_width * i / 2;  /* range / 2 */
        }

        for( j = i; j < class_count; j++ ) 
        {
            /* Count rising and falling slopes */
            assert( sum < ( RFC_COUNTS_LIMIT - rfm[ j-i, j ] - rfm[ j, j-i ] ) );
            sum += rfm[ j-i, j   ];
            sum += rfm[ j,   j-i ];
        }

        rp[i] = sum;
    }

    return true;
}


/**
 * @brief      Calculate the damage from a range pair histogram
 *
 * @param      ctx             The rainflow context
 * @param[in]  rp              The range pair counts, may be NULL (Mind that
 *                             full_inc describes 1 cycle!)
 * @param[in]  Sa              The buffer for amplitudes respective rp, may be
 *                             NULL
 * @param[out] damage          The buffer for cumulated damage
 * @param      rp_calc_method  The rp calculate method
 *                             (RFC_RP_DAMAGE_CALC_METHOD_*)
 *
 * @return     true on success
 */
bool RFC_damage_from_rp( const void *ctx, const rfc_counts_t *rp, const rfc_value_t *Sa, double *damage, rfc_rp_damage_method_e rp_calc_method )
{
    const unsigned    from = 0;
          unsigned    class_count;
          double      D;             /* Cumulative damage */
          int         i, j;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !damage )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !rp )
    {
        rp = rfc_ctx->rp;
    }

    if( !rp || !class_count )
    {
        return false;
    }

    /* Sa must be sorted in ascending order */
    for( i = 1; Sa && i < (int)class_count; i++ )
    {
        if( Sa[i] < Sa[i-1] )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
    }

    D = 0.0;

    /* Calculate the "Miner Consequent" approach */
    if( rp_calc_method == RFC_RP_DAMAGE_CALC_METHOD_CONSEQUENT )
    {
        rfc_wl_param_s      wl;
        double              q     = rfc_ctx->wl_q;   /* Fatigue strength depression exponent */
        double              Sd    = rfc_ctx->wl_sd;  /* Fatigue strength SD for the unimpaired(!) part */
        double              Nd    = rfc_ctx->wl_nd;  /* Fatigue strength cycle count ND for the unimpaired(!) part */
        double              Sj    = Sd;              /* Current fatigue strength, impacted part */
        double              D_inv = 0.0;             /* The inverse damage */
        bool                ok    = true;

        /* Omission not allowed here! */
        if( rfc_ctx->wl_omission > 0.0 )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }

        /* Backup WL parameters */
        RFC_wl_param_get( rfc_ctx, &wl );

        /* Remove fatigue strength temporarily */
        rfc_ctx->wl_sd = 0.0;
        rfc_ctx->wl_nd = DBL_MAX;

        for( j = (int)class_count - 1; j >= -1 && ok; j-- )
        {
            double D_j = 0.0;  /* Damage for partial histogram */
            double Sa_j;       /* New fatigue strength */
            double weight;     /* Weight for partial histogram */

            /* Get the new degraded fatigue strength in Sa_j */
            if( j >= 0 )
            {
                Sa_j = Sa ? Sa[j] : AMPLITUDE( rfc_ctx, j );
            }
            else
            {
                Sa_j = 0.0;
            }

            /* Forward until amplitude is below fatigue strength SD for the unimpaired part */
            if( Sa_j >= Sd && Sd > 0.0 ) continue;

            /**
             *  The following lines of code doesn't reflect the original formula as in [6] chapter 3.2.9, formula 3.2-65.
             *  This realization resolves his approach to use partial damages instead of particular cycles from the histogram.
             *  This procedure benefits from inserting an abstraction layer, which permits to use any representation of 
             *  a Woehler curve with a given fatigue strength limit providing a damage per cycle representation. 
             *  Note that the Miner Consequent approach, in contrast to Original/Elementary/Modified, is unsuitable when using 
             *  an inappropriate Woehler curve, to estimate so called "pseudo damages".
             */

            /* Weighted damage */
            weight = pow( Sj / Sd, q ) - pow( Sa_j / Sd, q );
            Sj     = Sa_j;

            if( weight <= 0.0 ) continue;

            /* Calculate damage for partial histogram */
            for( i = (int)class_count - 1; i > j; i-- )
            {
                double Sa_i = Sa ? Sa[i] : AMPLITUDE( rfc_ctx, i );
                double D_i;
                
                if( !damage_calc_amplitude( rfc_ctx, Sa_i, &D_i ) )
                {
                    ok = false;
                    break;
                }

                D_j += D_i * rp[i];
            }

            if( D_j > 0.0 )
            {
                D_inv += weight / D_j;
            }
        }

        /* Restore WL parameters */
        RFC_wl_param_set( rfc_ctx, &wl );

        if( !ok ) return false;

        /* Get the final damage value */
        D = 1.0 / D_inv;
    }
    else if( rp_calc_method == RFC_RP_DAMAGE_CALC_METHOD_ELEMENTAR )
    {
        rfc_wl_param_s  wl;
        bool ok;

        (void)RFC_wl_param_get( rfc_ctx, &wl );

        rfc_ctx->wl_sd = 0.0;
        rfc_ctx->wl_nd = DBL_MAX;
        rfc_ctx->wl_k2 = rfc_ctx->wl_k;
        rfc_ctx->wl_q2 = rfc_ctx->wl_q2;

#if RFC_DAMAGE_FAST
        {
            rfc_ctx->damage_lut_inapt++;
            ok = RFC_damage_from_rp( rfc_ctx, rp, Sa, damage, RFC_RP_DAMAGE_CALC_METHOD_DEFAULT );
            rfc_ctx->damage_lut_inapt--;
        }
#else /*!RFC_DAMAGE_FAST*/
        ok = RFC_damage_from_rp( rfc_ctx, rp, Sa, damage, RFC_RP_DAMAGE_CALC_TYPE_DEFAULT );
#endif /*RFC_DAMAGE_FAST*/

        (void)RFC_wl_param_set( rfc_ctx, &wl );

        return ok;
    }
    else if( rp_calc_method == RFC_RP_DAMAGE_CALC_METHOD_MODIFIED)
    {
        rfc_wl_param_s wl;
        bool ok;

        (void)RFC_wl_param_get( rfc_ctx, &wl );

        rfc_ctx->wl_sd = 0.0;
        rfc_ctx->wl_nd = DBL_MAX;

#if RFC_DAMAGE_FAST
        {
            rfc_ctx->damage_lut_inapt++;
            ok = RFC_damage_from_rp( rfc_ctx, rp, Sa, damage, RFC_RP_DAMAGE_CALC_METHOD_DEFAULT );
            rfc_ctx->damage_lut_inapt--;
        }
#else /*!RFC_DAMAGE_FAST*/
        ok = RFC_damage_from_rp( rfc_ctx, rp, Sa, damage, RFC_RP_DAMAGE_CALC_TYPE_DEFAULT );
#endif /*RFC_DAMAGE_FAST*/

        (void)RFC_wl_param_set( rfc_ctx, &wl );

        return ok;
    }
    else if( rp_calc_method == RFC_RP_DAMAGE_CALC_METHOD_DEFAULT )
    {
        for( i = 0; i < (int)class_count; i++ )
        {
            if( rp[i] )
            {
                double D_i;

                if( Sa )
                {
                    if( !damage_calc_amplitude( rfc_ctx, Sa[i], &D_i ) )
                    {
                        return false;
                    }
                    D += D_i * rp[i];
                }
                else
                {
                    if( !damage_calc( rfc_ctx, from, i /*to*/, &D_i, NULL /*Sa_ret*/ ) )
                    {
                        return false;
                    }
                    D += D_i * rp[i];
                }
            }
        }
    }
    else return false;

    *damage = D / rfc_ctx->full_inc;
    return true;
}


/**
 * @brief      Calculate the damage from rainflow matrix
 *
 * @param      ctx     The rainflow context
 * @param[in]  rfm     The input rainflow matrix (may be NULL), .full_inc represents one "cycle"!
 * @param[out] damage  The buffer for cumulated damage
 *
 * @return     true on success
 */
bool RFC_damage_from_rfm( const void *ctx, const rfc_counts_t *rfm, double *damage )
{
    unsigned    from, to;
    unsigned    class_count;
    double      D;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !damage )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    class_count = rfc_ctx->class_count;

    if( !rfm )
    {
        rfm = rfc_ctx->rfm;
    }

    if( !rfm || !class_count )
    {
        return false;
    }

    D = 0.0;
    for( from = 0; from < class_count; from++ )
    {
        for( to = 0; to < class_count; to++ )
        {
            if( rfm[ from * class_count + to ] )
            {
                double D_i;

                if( !damage_calc( rfc_ctx, from, to, &D_i, NULL /*Sa_ret*/ ) )
                {
                    return false;
                }
                D += D_i * rfm[ from * class_count + to ];
            }
        }
    }

    *damage = D / rfc_ctx->full_inc;
    return true;
}


/**
 * @brief      Calculate junction point between k and k2 for a Woehler curve
 *
 * @param      ctx   The rainflow context
 * @param      s0    Any point on slope k
 * @param      n0    Any point on slope k
 * @param      k     Slope k
 * @param[out] sx    Junction point between k and k2
 * @param      nx    Junction point between k and k2
 * @param      k2    Slope k2
 * @param      sd    The fatigue strength
 * @param      nd    The cycle count at sd
 *
 * @return     true on success
 */
bool RFC_wl_calc_sx( const void *ctx, double s0, double n0, double k, double *sx, double nx, double k2, double sd, double nd )
{
    double nom, den;

    k  = fabs(k);
    k2 = fabs(k2);

    if(    s0 <= 0.0    ||    n0 <= 0.0    ||
        /* sx <= 0.0    || */ nx <= 0.0    ||
           sd <= 0.0    ||    nd <= 0.0    ||    !sx )
    {
        return false;
    }

    nom = log(s0)*k - log(sd)*k2 + log(n0) - log(nd);
    den = k - k2;

    if( den == 0.0 ) return false;

    *sx = exp( nom / den );

    return true;
}


/**
 * @brief      Calculate fatigue strength for a Woehler curve
 *
 * @param      ctx   The rainflow context
 * @param      s0    Any point on slope k
 * @param      n0    Any point on slope k
 * @param      k     Slope k
 * @param      sx    Junction point between k and k2
 * @param      nx    The cycle
 * @param      k2    Slope k2
 * @param[out] sd    The fatigue strength
 * @param      nd    The cycle count at sd
 *
 * @return     true on success
 */
bool RFC_wl_calc_sd( const void *ctx, double s0, double n0, double k, double sx, double nx, double k2, double *sd, double nd )
{
    double nom, den;

    k  = fabs(k);
    k2 = fabs(k2);

    if(    s0 <= 0.0    ||    n0 <= 0.0    ||
           sx <= 0.0    ||    nx <= 0.0    ||
        /* sd <= 0.0    || */ nd <= 0.0    ||    !sd )
    {
        return false;
    }

    nom = log(s0)*k - log(sx)*(k-k2) + log(n0) - log(nd);
    den = k2;

    if( den == 0.0 ) return false;

    *sd = exp( nom / den );

    return true;
}


/**
 * @brief      Calculate the slope k2 for a Woehler curve
 *
 * @param      ctx   The rainflow context
 * @param      s0    Any point on slope k
 * @param      n0    Any point on slope k
 * @param      k     Slope k
 * @param      sx    Junction point between k and k2
 * @param      nx    Junction point between k and k2
 * @param[out] k2    Slope k2
 * @param      sd    The fatigue strength
 * @param      nd    The cycle count at sd
 *
 * @return     true on success
 */
bool RFC_wl_calc_k2( const void *ctx, double s0, double n0, double k, double sx, double nx, double *k2, double sd, double nd )
{
    double nom, den;

    k  = fabs(k);

    if(    s0 <= 0.0    ||    n0 <= 0.0    ||
           sx <= 0.0    ||    nx <= 0.0    ||
           sd <= 0.0    ||    nd <= 0.0    ||    !k2 )
    {
        return false;
    }

    nom = ( log(s0) - log(sx) ) * k + log(n0) - log(nd);
    den =   log(sd) - log(sx);

    if( den == 0.0 )
    {
        *k2 = -DBL_MAX;  /* Negative infinity */
    }
    else
    {
        *k2 = -fabs( nom / den );
    }

    return true;
}

/**
 * @brief      Calculate a point on a Woehler slope
 *
 * @param      ctx   The rainflow context
 * @param      s0    Any point on the slope, Sa
 * @param      n0    Any point on the slope, N
 * @param      k     The Woehler slope
 * @param      n     Given cycles for required Sa
 * @param[out] sa    The buffer for Sa(N)
 *
 * @return     true on success
 */
bool RFC_wl_calc_sa( const void *ctx, double s0, double n0, double k, double n, double *sa )
{
    /* (s0/Sa)^-k = n0/n */
    /* (Sa/s0)^k  = n0/n */

    k  = fabs(k);

    if(    s0 <= 0.0    ||    n0 <= 0.0    ||
           n  <= 0.0    ||   !sa )
    {
        return false;
    }

    *sa = pow( n0 / n, 1.0 / k ) * s0;

    return true;
}


/**
 * @brief      Calculate a point on a Woehler slope
 *
 * @param      ctx   The rainflow context
 * @param      s0    Any point on the slope, Sa
 * @param      n0    Any point on the slope, N
 * @param      k     The Woehler slope
 * @param      sa    Given amplitude Sa for required N
 * @param[out] n     The buffer for N(Sa)
 *
 * @return     true on success
 */
bool RFC_wl_calc_n( const void *ctx, double s0, double n0, double k, double sa, double *n )
{
    /* (Sa/s0)^-k = n/n0 */
    /* (s0/Sa)^k  = n/n0 */

    k  = fabs(k);

    if(    s0 <= 0.0    ||    n0 <= 0.0    ||
           sa <= 0.0    ||   !n )
    {
        return false;
    }

    *n = pow( s0 / sa, k ) * n0;

    return true;
}


/**
 * @brief      Set class parameter.
 *
 * @param      ctx          The rainflow context
 * @param[in]  class_param  The new class parameter
 *
 * @return     true on success
 */
bool RFC_class_param_set( void *ctx, const rfc_class_param_s *class_param )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }

    if( !class_param                          || 
         rfc_ctx->state     != RFC_STATE_INIT || 
         class_param->count <= 0              || 
         class_param->width <= 0.0            )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    rfc_ctx->class_count  = class_param->count;
    rfc_ctx->class_width  = class_param->width;
    rfc_ctx->class_offset = class_param->offset;

    return true;
}


/**
 * @brief      Get class parameter.
 *
 * @param      ctx          The rainflow context
 * @param[out] class_param  The class parameter
 *
 * @return     true on success
 */
bool RFC_class_param_get( const void *ctx, rfc_class_param_s *class_param )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( !class_param                     ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    class_param->count  = rfc_ctx->class_count;
    class_param->width  = rfc_ctx->class_width;
    class_param->offset = rfc_ctx->class_offset;

    return true;
}


/**
 * @brief      Get the class number for one value
 *
 * @param[in]  ctx     The rainflow context
 * @param[in]  value   The value
 * @param      number  The class number (base 0)
 *
 * @return     true on success
 */
bool RFC_class_number( const void *ctx, rfc_value_t value, unsigned *class_number )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( !class_number                    ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    *class_number = QUANTIZE( rfc_ctx, value );

    return true;
}


/**
 * @brief      Get the class upper value for one class
 *
 * @param[in]  ctx                The rainflow context
 * @param[in]  class_number       The class number (base 0)
 * @param      class_upper_value  The class upper value
 *
 * @return     true on success
 */
bool RFC_class_upper( const void *ctx, unsigned class_number, rfc_value_t *class_upper_value )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( !class_upper_value               ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    *class_upper_value = CLASS_UPPER( rfc_ctx, class_number );

    return true;
}


/**
 * @brief      Get class mean value for one class
 *
 * @param[in]  ctx               The rainflow context
 * @param[in]  class_number      The class number
 * @param      class_mean_value  The class mean value
 *
 * @return     true on success
 */
bool RFC_class_mean( const void *ctx, unsigned class_number, rfc_value_t *class_mean_value )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( !class_mean_value                ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    *class_mean_value = CLASS_MEAN( rfc_ctx, class_number );

    return true;
}


/**
 * @brief      Set flags
 *
 * @param      ctx    The rainflow context
 * @param[in]  flags  The flags
 * @param[in]  stack  ID of flags stack
 *
 * @return     true on success
 */
bool RFC_set_flags( void *ctx, int flags, int stack )
{
    RFC_CTX_CHECK_AND_ASSIGN

    switch( stack )
    {
        case 0:
            rfc_ctx->internal.flags = (rfc_flags_e)flags;
            break;

#if RFC_DEBUG_FLAGS
        case 1:
            rfc_ctx->internal.debug_flags = (rfc_debug_flags_e)flags;
            break;
#endif /*RFC_DEBUG_FLAGS*/

        default:
            return false;
    }

    return true;
}


/**
 * @brief      Get flags
 *
 * @param      ctx    The rainflow context
 * @param[out] flags  Flag for debugging flags
 * @param[in]  stack  ID of flags stack
 *
 * @return     true on success
 */
bool RFC_get_flags( const void *ctx, int *flags, int stack )
{
    RFC_CTX_CHECK_AND_ASSIGN

    switch( stack )
    {
        case 0:
            *flags = (int)rfc_ctx->internal.flags;
            break;

#if RFC_DEBUG_FLAGS
        case 1:
            *flags = (int)rfc_ctx->internal.debug_flags;
            break;
#endif /*RFC_DEBUG_FLAGS*/

        default:
            return false;
    }

    return true;
}


/**
 * @brief      Set Woehler curve parameters
 *
 * @param      ctx       The rainflow context
 * @param[in]  wl_param  The Woehler curve parameters
 *
 * @return     true on success
 */
bool RFC_wl_param_set( void *ctx, const rfc_wl_param_s *wl_param )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }
    
    if( !wl_param                        ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    rfc_ctx->wl_sx          = wl_param->sx;
    rfc_ctx->wl_nx          = wl_param->nx;
    rfc_ctx->wl_k           = wl_param->k;
    rfc_ctx->wl_q           = wl_param->q;
    rfc_ctx->wl_sd          = wl_param->sd;
    rfc_ctx->wl_nd          = wl_param->nd;
    rfc_ctx->wl_k2          = wl_param->k2;
    rfc_ctx->wl_q2          = wl_param->q2;
    rfc_ctx->wl_omission    = wl_param->omission;

    return true;
}


/**
 * @brief      Get Woehler curve parameters
 *
 * @param      ctx       The rainflow context
 * @param[out] wl_param  The Woehler curve parameters
 *
 * @return     true on success
 */
bool RFC_wl_param_get( const void *ctx, rfc_wl_param_s *wl_param )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state > RFC_STATE_FINISHED )
    {
        return false;
    }
    
    if( !wl_param                        ||
         rfc_ctx->state < RFC_STATE_INIT )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }

    wl_param->sx            = rfc_ctx->wl_sx;
    wl_param->nx            = rfc_ctx->wl_nx;
    wl_param->k             = rfc_ctx->wl_k;
    wl_param->q             = rfc_ctx->wl_q;
    wl_param->sd            = rfc_ctx->wl_sd;
    wl_param->nd            = rfc_ctx->wl_nd;
    wl_param->q2            = rfc_ctx->wl_q2;
    wl_param->k2            = rfc_ctx->wl_k2;
    wl_param->omission      = rfc_ctx->wl_omission;
    wl_param->D             = 0.0;

    return true;
}
#endif /*!RFC_MINIMAL*/



#if RFC_AT_SUPPORT
/**
 * @brief      Amplitude transformation to take mean load influence into
 *             account.
 *
 * @param      ctx             The rainflow context
 * @param      Sa              Amplitude
 * @param      Sm              Mean load
 * @param[out] Sa_transformed  Transformed amplitude Sa
 *
 * @return     true on success
 */
bool RFC_at_transform( const void *ctx, double Sa, double Sm, double *Sa_transformed )
{
    double Sa_transform = Sa;

    RFC_CTX_CHECK_AND_ASSIGN

    if( !Sa_transformed )
    {
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );
    }
    
    /* Amplitude is always positive */
    Sa = fabs( Sa );

#if RFC_USE_DELEGATES
    if( rfc_ctx->at_transform_fcn )
    {
        return rfc_ctx->at_transform_fcn( rfc_ctx, Sa, Sm, Sa_transformed );
    }
#endif

    if( !rfc_ctx->at.count )
    {
        /* No reference curve given, return original amplitude */
        *Sa_transformed = Sa;
        return true;
    }

    if( Sa == 0.0 )
    {
        /* Zero amplitude */
        Sa_transform = 0.0;
    }
    else
    {
        double Sm_norm_base;
        double Sm_norm_target;
        double alleviation_base;
        double alleviation_target;

        /* Normalize Sm (Sa=1) */
        Sm_norm_base = Sm / Sa;
        if( !at_alleviation( rfc_ctx, Sm_norm_base, &alleviation_base ) )
        {
            return false;
        }

        if( rfc_ctx->at.R_pinned )
        {
            /* Calculate intersection of R slope and M slope */
            if( !at_R_to_Sm_norm( rfc_ctx, rfc_ctx->at.R_rig, &Sm_norm_target ) )
            {
                return false;
            }

            if( !at_alleviation( rfc_ctx, Sm_norm_target, &alleviation_target ) )
            {
                return false;
            }

            Sa_transform = Sa / alleviation_base * alleviation_target;
        }
        else
        {
            /* Calculate intersection of mean load on test rig (Sm_rig) with curve taken from Haigh diagram */
            const double   *Sa_   = rfc_ctx->at.Sa;
            const double   *Sm_   = rfc_ctx->at.Sm;
                  unsigned  count = rfc_ctx->at.count;
                  unsigned  n;
                  double    Sa_lhs, Sa_rhs;
                  double    Sm_lhs, Sm_rhs;

            /* Check each segment from the reference curve */
            for( n = 0; n <= count; n++ )
            {
                if( n )
                {
                    Sa_lhs = Sa_rhs;
                    Sm_lhs = Sm_rhs;

                    if( n < count )
                    {
                        if( !at_alleviation( rfc_ctx, Sm_[n] / Sa_[n], &alleviation_target ) )
                        {
                            return false;
                        }

                        /* Next segment */
                        Sa_rhs = Sa / alleviation_base * alleviation_target;
                        Sm_rhs = Sa_rhs / Sa_[n] * Sm_[n];
                    }
                    else
                    {
                        /* Last segment */
                        assert( Sm_lhs <= rfc_ctx->at.Sm_rig );

                        Sa_transform = Sa_lhs;
                    }
                }
                else /* n == 0 */
                {
                    if( !at_alleviation( rfc_ctx, Sm_[0] / Sa_[0], &alleviation_target ) )
                    {
                        return false;
                    }

                    /* First segment */
                    Sa_rhs = Sa / alleviation_base * alleviation_target;
                    Sm_rhs = Sa_rhs / Sa_[0] * Sm_[0];

                    if( rfc_ctx->at.Sm_rig <= Sm_rhs)
                    {
                        Sa_transform = Sa_rhs;
                        break;
                    }
                    else continue;
                }

                /* Interpolate intermediate points */
                if( Sm_lhs <= rfc_ctx->at.Sm_rig && rfc_ctx->at.Sm_rig <= Sm_rhs )
                {
                    double frac, denom;

                    denom = Sm_rhs - Sm_lhs;

                    frac = ( denom < 1e-20 ) ? 1.0 : ( ( rfc_ctx->at.Sm_rig - Sm_lhs ) / denom );

                    Sa_transform = Sa_lhs * ( 1.0 - frac ) + Sa_rhs * frac;
                    break;
                }
            }
        }
    }

    *Sa_transformed = Sa_transform;

    return true;
}
#endif /*RFC_AT_SUPPORT*/

#if RFC_DEBUG_FLAGS
int RFC_debug_fprintf( void *ctx, FILE *stream, const char *fmt, ... )
{
    int result;
    va_list arg;

    RFC_CTX_CHECK_AND_ASSIGN

    va_start( arg, fmt );

#if RFC_USE_DELEGATES
    if( rfc_ctx->debug_vfprintf_fcn )
    {
        result = rfc_ctx->debug_vfprintf_fcn( rfc_ctx, stream, fmt, arg );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        result = vfprintf( stream, fmt, arg );
    }

    va_end( arg );

    return result;
}
#endif /*RFC_DEBUG_FLAGS*/










/*** Implementation static functions ***/

#if RFC_AT_SUPPORT
/**
 * @brief      Clear all data for amplitude transformation
 *
 * @param      rfc_ctx  The rainflow context
 */
static
void clear_at( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx );

    rfc_ctx->at.Sa                      = NULL;
    rfc_ctx->at.Sm                      = NULL;
    rfc_ctx->at.count                   = 0;
    rfc_ctx->at.M                       = 0.0;
    rfc_ctx->at.Sm_rig                  = 0.0;
    rfc_ctx->at.R_rig                   = 0.0;
    rfc_ctx->at.R_pinned                = false;

    rfc_ctx->internal.at_haigh.count    = 0;

}
#endif /*RFC_AT_SUPPORT*/


#if RFC_DAMAGE_FAST
static
void clear_lut( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state >= RFC_STATE_INIT );

    if( rfc_ctx->damage_lut )
    {
        memset( rfc_ctx->damage_lut, 0, sizeof(double) * rfc_ctx->class_count * rfc_ctx->class_count );
    }
    rfc_ctx->damage_lut_inapt = 1;

#if RFC_AT_SUPPORT
    if( rfc_ctx->amplitude_lut )
    {
        memset( rfc_ctx->amplitude_lut, 0, sizeof(double) * rfc_ctx->class_count * rfc_ctx->class_count );
    }
#endif /*RFC_AT_SUPPORT*/
}
#endif /*RFC_DH_SUPPORT*/


/**
 * @brief      Processing one data point. Find turning points and check for
 *             closed cycles.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  pt       The data tuple
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool feed_once( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt, rfc_flags_e flags )
{
    rfc_value_tuple_s *tp_residue;  /* Pointer to residue element */

    assert( rfc_ctx && pt );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

#if RFC_DH_SUPPORT
    /* Resize damage history if necessary */
    if( !feed_once_dh( rfc_ctx, pt ) )
    {
        return false;
    }
#endif /*RFC_DH_SUPPORT*/

    /* Check for next turning point and update residue. tp_residue is NULL, if there is no turning point */
    /* Otherwise tp_residue refers the forelast element in member rfc_ctx->residue */
    tp_residue = feed_filter_pt( rfc_ctx, pt );

#if RFC_TP_SUPPORT
    /* Check if pt influences margins (tp_residue may be set to NULL then!) */
    if( !feed_once_tp_check_margin( rfc_ctx, pt, &tp_residue ) )
    {
        return false;
    }
#endif /*RFC_TP_SUPPORT*/

    /* Countings */

    /* Add turning point and check for closed cycles */
    if( tp_residue )
    {
#if RFC_TP_SUPPORT
        /* Add a copy of tp_residue to rfc_ctx->tp and alter tp_residue->tp_pos to its position in rfc_ctx->tp */
        if( !tp_set( rfc_ctx, 0, tp_residue ) )
        {
            return false;
        }
#endif /*RFC_TP_SUPPORT*/

#if !RFC_MINIMAL
        /* New turning point, do LC count */
        cycle_process_lc( rfc_ctx, flags & (RFC_FLAGS_COUNT_LC | RFC_FLAGS_ENFORCE_MARGIN) );
        flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

        if( rfc_ctx->class_count )
        {
            /* Check for closed cycles and count. Modifies residue! */
            cycle_find( rfc_ctx, flags );
        }
        else
        {
            if( rfc_ctx->residue_cnt > 1 )
            {
                residue_remove_item( rfc_ctx, 0, 1 );
            }
        }
    }

    return true;
}


#if RFC_DH_SUPPORT
/**
 * @brief      Resize damage history if necessary.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  pt       The point
 *
 * @return     true on success
 */
static
bool feed_once_dh( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt )
{
    assert( rfc_ctx );
    assert( pt );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    if( rfc_ctx->dh )
    {
        if( pt->pos > rfc_ctx->dh_cap )
        {
            size_t new_cap = (size_t)1024 * ( pt->pos / 640 + 1 ); /* + 60% + 1024 */

            rfc_ctx->dh = (double*)rfc_ctx->mem_alloc( rfc_ctx->dh, new_cap, 
                                                       sizeof(rfc_value_t), RFC_MEM_AIM_DH );

            if( !rfc_ctx->dh )
            {
                return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
            }

            memset( rfc_ctx->dh + rfc_ctx->dh_cnt, 0, sizeof(rfc_value_t) * ( new_cap - rfc_ctx->dh_cap ) );
            rfc_ctx->dh_cap = new_cap;
        }

        rfc_ctx->dh_cnt = pt->pos;
    }

    return true;
}
#endif /*RFC_DH_SUPPORT*/


#if RFC_TP_SUPPORT
/**
 * @brief         Check if pt influences margins.
 *
 * @param         rfc_ctx     The rainflow context
 * @param[in]     pt          The new data point
 * @param[in,out] tp_residue  The new turning point (or NULL)
 *
 * @return        true on success
 */
bool feed_once_tp_check_margin( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt, rfc_value_tuple_s** tp_residue )
{
    bool do_margin;

    assert( rfc_ctx && tp_residue );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    do_margin = rfc_ctx->internal.flags & RFC_FLAGS_ENFORCE_MARGIN;

    if( do_margin && !rfc_ctx->tp_locked )
    {
        /* Check for left and right margin */
        switch( rfc_ctx->internal.margin_stage )
        {
            case 0:
            {
                rfc_value_tuple_s pt_left = *pt;

                assert( *tp_residue == NULL );

                /* Save left margin */
                rfc_ctx->internal.margin[0]  = *pt;

                /* Enqueue left margin as turning point */
                if( !tp_set( rfc_ctx, 0, &pt_left ) ) return false;

                rfc_ctx->internal.margin_stage = 1;
                break;
            }

            case 1:
            {
                /* Save right margin so far */
                rfc_ctx->internal.margin[1] = *pt;

                /* Prevent storing 1st point twice */
                if( *tp_residue )
                {
                    /* First turning point found */
                    rfc_ctx->internal.margin_stage = 2;

                    if( (*tp_residue)->value == rfc_ctx->internal.margin[0].value )
                    {
                        assert( rfc_ctx->tp_cnt == 1 );

                        /* Left margin and first turning point are identical, set reference in residue */
                        (*tp_residue)->tp_pos = 1;
                        (*tp_residue) = 0;  /* Avoid further processing of this turning point */
                    }
                }
                break;
            }

            case 2:
                /* Save right margin so far */
                rfc_ctx->internal.margin[1] = *pt;
                break;

            default:
                assert( false );
        }
    }

    return true;
}
#endif /*RFC_TP_SUPPORT*/


/**
 * @brief      Handling interim turning point and margin. If there are still
 *             unhandled turning point left, "finalizing" takes this into
 *             account for the rainflow algorithm.
 *
 * @param      rfc_ctx  The rainflow context
 *
 * @return     true on success
 */
static
bool feed_finalize( rfc_ctx_s *rfc_ctx )
{
    rfc_value_tuple_s *tp_interim = NULL;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    if( rfc_ctx->state < RFC_STATE_FINALIZE )
    {
        /* Adjust residue: Incorporate interim turning point */
        if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
        {
            tp_interim = &rfc_ctx->residue[rfc_ctx->residue_cnt];
            rfc_ctx->residue_cnt++;

            rfc_ctx->state = RFC_STATE_BUSY;
        }

#if RFC_TP_SUPPORT
        /* Finalize turning point storage */
        if( !feed_finalize_tp( rfc_ctx, tp_interim, rfc_ctx->internal.flags ) )
        {
            return false;
        }
#endif /*RFC_TP_SUPPORT*/

        if( tp_interim )
        {
            int flags = rfc_ctx->internal.flags;
#if !RFC_MINIMAL
            /* New turning point, do LC count */
            cycle_process_lc( rfc_ctx, rfc_ctx->internal.flags & (RFC_FLAGS_COUNT_LC | RFC_FLAGS_ENFORCE_MARGIN) );
            flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

            /* Check once more if a new cycle is closed now */
            cycle_find( rfc_ctx, flags );
        }

#if RFC_HCM_SUPPORT
        if( !feed_finalize_hcm( rfc_ctx, rfc_ctx->internal.flags ) )
        {
            return false;
        }
#endif /*RFC_HCM_SUPPORT*/

        rfc_ctx->state = RFC_STATE_FINALIZE;
    }

    return true;
}


#if RFC_TP_SUPPORT
/**
 * @brief         Finalize turning point storage.
 *
 * @param         rfc_ctx     The rainflow context
 * @param[in,out] tp_interim  The interim turning point (or NULL)
 * @param         flags       Only flag RFC_FLAGS_ENFORCE_MARGIN encounters
 *
 * @return        true on success
 */
static
bool feed_finalize_tp( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *tp_interim, rfc_flags_e flags )
{
    bool do_margin;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );
    
    /* Finalize turning points storage */
    do_margin = rfc_ctx->internal.flags & RFC_FLAGS_ENFORCE_MARGIN;
    if( do_margin && !rfc_ctx->tp_locked )
    {
        rfc_value_tuple_s *pt_right = &rfc_ctx->internal.margin[1];

        if( tp_interim )
        {
            if( rfc_ctx->internal.margin_stage > 0 && tp_interim->value == pt_right->value )
            {
                if( !tp_set( rfc_ctx, 0, pt_right ) ) return false;
                tp_interim->tp_pos = pt_right->tp_pos;
            }
            else
            {
                if( !tp_set( rfc_ctx, 0, tp_interim ) ) return false;
                if( !tp_set( rfc_ctx, 0, pt_right ) )   return false;
            }
        }
        else if( pt_right->pos > 0 )
        {
            if( !tp_set( rfc_ctx, 0, pt_right ) )   return false;
        }
    }
    else if( tp_interim )
    {
        if( !tp_set( rfc_ctx, 0, tp_interim ) ) return false;
    }

    /* Lock turning points queue */
    tp_lock( rfc_ctx, true );
    return true;
}
#endif /*RFC_TP_SUPPORT*/


#if RFC_HCM_SUPPORT
/**
 * @brief      Finalize HCM algorithm, copy residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool feed_finalize_hcm( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Move HCM stack to residue */
    if( rfc_ctx->counting_method == RFC_COUNTING_METHOD_HCM )
    {
        int stack_cnt = rfc_ctx->internal.hcm.IZ; /* Number of turning points in HCM stack */

        if( stack_cnt )
        {
            /* Reallocate residue */
            rfc_ctx->residue = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( rfc_ctx->residue, (size_t)stack_cnt, 
                                                                       sizeof(rfc_value_tuple_s), RFC_MEM_AIM_RESIDUE );

            if( !rfc_ctx->residue )
            {
                return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
            }

            memcpy( rfc_ctx->residue, rfc_ctx->internal.hcm.stack, sizeof(rfc_value_tuple_s) * stack_cnt );

            rfc_ctx->residue_cap = stack_cnt;
            rfc_ctx->residue_cnt = stack_cnt;

            /* Make HCM stack empty */
            rfc_ctx->internal.hcm.IZ = 0;
            rfc_ctx->internal.hcm.IR = 1;
        }
    }

    return true;
}
#endif /*RFC_HCM_SUPPORT*/


/**
 * @brief      Finalize pending counts, ignore residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_ignore( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Include interim turning point */
    return feed_finalize( rfc_ctx );
}


#if !RFC_MINIMAL
/**
 * @brief      Finalize pending counts, discard residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_discard( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Include interim turning point */
    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Finalize pending counts, weight unclosed cycles into RP,
 *             RF-matrix and damage and discard residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      weight   The weight for closed cycles
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_weight_cycles( rfc_ctx_s *rfc_ctx, rfc_counts_t weight, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Include interim turning point */
    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* Count every unclosed cycle with the given weight */
    if( rfc_ctx->residue && rfc_ctx->residue_cnt >= 2 )
    {
        size_t             i;
        int                flags    = rfc_ctx->internal.flags;
        rfc_value_tuple_s *from     = rfc_ctx->residue;
        rfc_counts_t       old_inc  = rfc_ctx->curr_inc;

        rfc_ctx->curr_inc = weight;

        for( i = 0; i + 1 < rfc_ctx->residue_cnt; i++ )
        {
            rfc_value_tuple_s *to   = from + 1;
            rfc_value_tuple_s *next = ( i + 2 < rfc_ctx->residue_cnt ) ? to + 1 : NULL;

            cycle_process_counts( rfc_ctx, from, to, next, flags );

            from = to;
        }

        rfc_ctx->curr_inc = old_inc;
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Finalize pending counts to fit HCM results.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_clormann_seeger( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    size_t i;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Include interim turning point */
    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    if( rfc_ctx->counting_method == RFC_COUNTING_METHOD_4PTM )
    {
        /* Counting correction on residue */

        for( i = 0; i + 4 < rfc_ctx->residue_cnt; )
        {
            size_t idx = rfc_ctx->residue_cnt + i;

            double A = (double)rfc_ctx->residue[idx+0].value;
            double B = (double)rfc_ctx->residue[idx+1].value;
            double C = (double)rfc_ctx->residue[idx+2].value;
            double D = (double)rfc_ctx->residue[idx+3].value;

            if( B * C < 0.0 && fabs(D) >= fabs(B) && fabs(B) >= fabs(C) )
            {
                rfc_value_tuple_s *from = &rfc_ctx->residue[idx+1];
                rfc_value_tuple_s *to   = &rfc_ctx->residue[idx+2];

                cycle_process_counts( rfc_ctx, from, to, to + 1, flags );

                /* Remove two inner turning points (idx+1 and idx+2) */
                residue_remove_item( rfc_ctx, i + 1, 2 );
                rfc_ctx->residue_cnt -= 2;
            }
            else
            {
                i++;
            }
        }
    }

    return true;
}


/**
 * @brief      Finalize pending counts, DIN method.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_rp_DIN45667( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Include interim turning point */
    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    if( flags && rfc_ctx->residue_cnt > 2 )
    {
        int i, j, k;
        int slopes_cnt = (int)rfc_ctx->residue_cnt - 1;
        struct slopes
        {
            int slope;
            rfc_value_tuple_s *lhs;
            rfc_value_tuple_s *rhs;
        } *slopes, tmp;

        slopes = (struct slopes*)rfc_ctx->mem_alloc( NULL, slopes_cnt, sizeof( slopes[0] ), RFC_MEM_AIM_TEMP );

        if( !slopes )
        {
            return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }

        /* Evaluate slopes */
        for( i = 0; i < slopes_cnt; i++ )
        {
            slopes[i].lhs   = &rfc_ctx->residue[i];
            slopes[i].rhs   = &rfc_ctx->residue[i+1];
            slopes[i].slope = (int)slopes[i].rhs->cls - slopes[i].lhs->cls;
        }

        /* Order slopes (simple "bubble sort") */
        k = 0;
        for( i = 0; i < slopes_cnt; i++ )
        {
            /* Drag higher ranks */
            for( j = slopes_cnt - 2; j >= i; j-- )
            {
                bool do_swap = false;

                if( ( slopes[j].slope > 0 ) == ( slopes[j+1].slope > 0 ) )
                {
                    do_swap = abs( slopes[j].slope ) < abs( slopes[j+1].slope );
                }
                else
                {
                    do_swap = slopes[j].slope < 0;
                    k       = j + 1;  /* k indicates the first falling slope now */
                }

                if( do_swap )
                {
                    tmp         = slopes[j];
                    slopes[j]   = slopes[j+1];
                    slopes[j+1] = tmp;
                }
            }
        }

        assert( k > 0 );

        /* Compare positive slopes with adjacent negative slopes */
        for( i = 0; i < k && i + k < slopes_cnt; i++ )
        {
            /* Get the smaller slope of two adjacent. On equality choose the primarily in time */
            int lh_slope = slopes[i].slope;
            int rh_slope = slopes[i+k].slope;

            if( abs( lh_slope ) == abs( rh_slope ) )
            {
                j = ( slopes[i].rhs->pos < slopes[i+k].rhs->pos ) ? i : (i+k);
            }
            else
            {
                j = ( abs( lh_slope ) < abs( rh_slope ) ) ? i : (i+k);
            }

            /* Do the countings for the matching slope */
            cycle_process_counts( rfc_ctx, slopes[j].lhs, slopes[j].rhs, NULL, flags );
        }

        rfc_ctx->mem_alloc( slopes, 0, 0, RFC_MEM_AIM_TEMP );
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Finalize pending counts, repeated residue method.
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_repeated( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
#if RFC_TP_SUPPORT
    size_t tp_cnt;
#endif /*RFC_TP_SUPPORT*/

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

#if RFC_TP_SUPPORT
    /* Keep number of turning points. The only scenario in which the number increases, is where
       the interim turning point gets a real turning point */
    tp_cnt = rfc_ctx->tp_cnt;
#endif /*RFC_TP_SUPPORT*/

    if( rfc_ctx->residue && rfc_ctx->residue_cnt && flags )
    {
        /* Feed again with the content of the residue itself */
        size_t              cnt     = rfc_ctx->residue_cnt;
        rfc_value_tuple_s  *residue = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, cnt + 1, 
                                                                              sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TEMP );

        if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
        {
            cnt++;
        }

        if( residue )
        {
            /* Make a copy of the residue */
            bool                     ok;
            size_t                   n         = cnt; 
            int                      old_flags = rfc_ctx->internal.flags;
            const rfc_value_tuple_s *from      = rfc_ctx->residue;
                  rfc_value_tuple_s *to        = residue;

            while( n-- )
            {
                *to++ = *from++;
            }

            rfc_ctx->internal.flags = flags;
            /* Feed again with the copy, no new turning points are generated, since residue[].tp_pos */

#if RFC_TP_SUPPORT
            do
            {
                size_t tp_cnt = rfc_ctx->tp_cnt;

                ok = RFC_feed_tuple( rfc_ctx, residue, cnt - 1 );

                if( ok )
                {
                    if( rfc_ctx->tp_cnt > tp_cnt )
                    {
                        /* Interim turning point became a turning point */
                        residue[cnt-1].tp_pos = tp_cnt + 1;
                    }
                    ok = RFC_feed_tuple( rfc_ctx, residue + cnt - 1, 1 );
                }
            } while(0);
#else /*!RFC_TP_SUPPORT*/
            ok = RFC_feed_tuple( rfc_ctx, residue, cnt );
#endif /*RFC_TP_SUPPORT*/

            rfc_ctx->internal.flags = old_flags;

            /* Free temporary residue */
            rfc_ctx->mem_alloc( residue, 0, 0, RFC_MEM_AIM_TEMP );

            if( !ok )
            {
                return false;
            }
        }
        else
        {
            return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }
    }

    /* Include interim turning point */
    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief         Backup/restore of residue
 *
 * @param         rfc_ctx      The rainflow context
 * @param[in,out] residue      The copy of the current residue
 * @param[out]    residue_cap  The capacity of the given residue
 * @param[out]    residue_cnt  The number of points in the given residue
 * @param         restore      true->restore, false->backup
 *
 * @return        true on success
 */
static
bool residue_exchange( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s **residue, size_t *residue_cap, size_t *residue_cnt, bool restore )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );
    assert( residue && residue_cap && residue_cnt );
    assert( rfc_ctx->residue_cap > 0 );

    if( !restore )
    {
        /* Backup */
        *residue_cap = rfc_ctx->residue_cap;
        *residue_cnt = rfc_ctx->residue_cnt;
        *residue     = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, *residue_cap, 
                                                               sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TEMP );

        if( !*residue )
        {
            return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }

        memcpy( *residue, rfc_ctx->residue, *residue_cnt * sizeof(rfc_value_tuple_s) );
    }
    else
    {
        /* Restore */

        /* Release residue */
        mem_alloc( rfc_ctx->residue, /*num*/ 0, /*size*/ 0, RFC_MEM_AIM_TEMP );

        /* Assign backup */
        rfc_ctx->residue_cap = *residue_cap;
        rfc_ctx->residue_cnt = *residue_cnt;
        rfc_ctx->residue     = *residue;
    }

    return true;
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Remove items (points) from the residue
 *
 * @param      rfc_ctx  The rainflow context
 * @param      index    The item position in residue, base 0
 * @param      count    The number of elements to remove
 */
static
void residue_remove_item( rfc_ctx_s *rfc_ctx, size_t index, size_t count )
{
    size_t  from = index + count,
            to   = index, 
            end;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );
    assert( rfc_ctx->residue && index + count <= rfc_ctx->residue_cnt );

    end = (int)rfc_ctx->residue_cnt;

    /* Example
                         |cnt(5)
                         v
               |O|O|X|O|O|
                   ^
                   |index(2)

        elements = 2    (5-2-1+0 = 2, no interim turning point)
    */

    if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
    {
        /* Include interim turning point */
        end++;
    }

    /* Shift points */
    while( from < end )
    {
        rfc_ctx->residue[to++] = rfc_ctx->residue[from++];
    }

    rfc_ctx->residue_cnt -= count;
}


/**
 * @brief      Calculate damage for one cycle with given amplitude Sa
 *
 * @param      rfc_ctx  The rainflow context
 * @param      Sa       The amplitude
 * @param[out] damage   The damage
 *
 * @return     true on success
 */
static
bool damage_calc_amplitude( rfc_ctx_s *rfc_ctx, double Sa, double *damage )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    do {
        /* Damage */
        double D = 0.0;

        if( Sa >= 0.0 )
        {
            /* D =           h /    ND   *    ( Sa /    SD)  ^ ABS(k)   */
            /* D = exp(  log(h /    ND)  + log( Sa /    SD)  * ABS(k) ) */
            /* D = exp( (log(h)-log(ND)) + (log(Sa)-log(SD)) * ABS(k) ) */
            /* D = exp(      0 -log(ND)  + (log(Sa)-log(SD)) * ABS(k) ) */

#if !RFC_MINIMAL
            if( Sa > rfc_ctx->wl_omission )
            {
                /* Constants for the Woehler curve */
                const double SX_log = log(rfc_ctx->wl_sx);
                const double NX_log = log(rfc_ctx->wl_nx);

                if( Sa > rfc_ctx->wl_sx )
                {
                    const double k = rfc_ctx->wl_k;
                    
                    /* Upper slope, finite life scope */
                    D = exp( fabs(k) * ( log(Sa) - SX_log ) - NX_log );
                }
                else if( Sa > rfc_ctx->wl_sd )
                {
                    const double k2 = rfc_ctx->wl_k2;

                    /* Lower slope, transition scope, modified Miners' rule */
                    D = exp( fabs(k2) * ( log(Sa) - SX_log ) - NX_log );
                }
                /* Amplitudes below fatigue strength have no influence */
            }
#else /*RFC_MINIMAL*/
            /* Constants for the Woehler curve */
            const double SX_log = log(rfc_ctx->wl_sx);
            const double NX_log = log(rfc_ctx->wl_nx);
            const double k      = rfc_ctx->wl_k;

            /* Miner original */
            D = exp( fabs(k)  * ( log(Sa) - SX_log ) - NX_log );
#endif /*!RFC_MINIMAL*/
        }
        else
        {
            assert( false );
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }

        *damage = D;
        
    } while(0);

    return true;
}


#if RFC_AT_SUPPORT
/**
 * @brief      Calculate the normalized mean load (Sa=1) for a given load ratio
 *             "R"
 *
 * @param      rfc_ctx  The rainflow context
 * @param      R        Load ratio (Su/So)
 * @param[out] Sm_norm  Normalized mean load
 *
 * @return     true on success
 */
static
bool at_R_to_Sm_norm( rfc_ctx_s *rfc_ctx, double R, double *Sm_norm )
{
    if( isinf( R ) || !Sm_norm )
    {
        return error_raise( rfc_ctx, RFC_ERROR_AT );
    }

    /* Su      = Sm - Sa */
    /* So      = Sm + Sa */
    /* R       = Su / So */
    /* Sm      = Sa * ( 1 + R ) / ( 1 - R ) */
    /* Sm_norm = Sm / Sa */

    *Sm_norm = ( 1 + R ) / ( 1 - R );

    return true;
}


/**
 * @brief      Calculate the influence of (normalized, Sa=1) mean load on
 *             fatigue strength (Haigh-diagram)
 *
 * @param      rfc_ctx      The rainflow context
 * @param      Sm_norm      Mean load, normalized (Sm/Sa)
 * @param[out] alleviation  Alleviation factor on reference curve
 *
 * @return     true on success
 */
static
bool at_alleviation( rfc_ctx_s *rfc_ctx, double Sm_norm, double *alleviation )
{
    assert( rfc_ctx );

    if( !alleviation )
    {
        return error_raise( rfc_ctx, RFC_ERROR_AT );
    }

    if( !rfc_ctx->at.count )
    {
        /* No reference curve given, no transformation */
        *alleviation = 1.0;
    }
    else
    {
        /* Reference curve */
        const double   *Sa_   = rfc_ctx->at.Sa;
        const double   *Sm_   = rfc_ctx->at.Sm;
              unsigned  count = rfc_ctx->at.count;

        assert( rfc_ctx && Sm_ && Sa_ );

        if( Sm_norm <= Sm_[0] / Sa_[0] )
        {
            /* First segment */
            *alleviation = Sa_[0];  /* Clip to first point */
        }
        else if( Sm_norm >= Sm_[count-1] / Sa_[count-1] )
        {
            /* Last segment */
            *alleviation = Sa_[count-1];  /* Clip to last point */
        }
        else
        {
            unsigned i;
            bool done = false;

            /* Select correct segment */
            for( i = 1; !done && i < count; i++ )
            {
                assert( Sa_[i-1] > 0.0 && Sa_[i] > 0.0 && Sm_[i-1] <= Sm_[i] );

                if( Sm_norm > Sm_[i-1] / Sa_[i-1] && Sm_norm <= Sm_[i] / Sa_[i] )
                {
                    /* Intersection of R slope and M slope */
                    double M_signed = ( Sa_[i] - Sa_[i-1] ) / ( Sm_[i] - Sm_[i-1] );

                    *alleviation = ( Sa_[i-1] - M_signed * Sm_[i-1] ) / ( 1.0 - M_signed * Sm_norm );
                    done = true;
                }
            }

            if( !done )
            {
                assert( false );
                return error_raise( rfc_ctx, RFC_ERROR_AT );
            }
        }
    }

    return true;
}
#endif /*RFC_AT_SUPPORT*/


/**
 * @brief      Calculate fictive damage for one closed (full) cycle.
 *
 * @param      rfc_ctx     The rainflow context
 * @param      class_from  The starting class
 * @param      class_to    The ending class
 * @param[out] damage      The damage value for the closed cycle
 * @param[out] Sa_ret      The amplitude, may be NULL
 *
 * @return     true on success
 */
static
bool damage_calc( rfc_ctx_s *rfc_ctx, unsigned class_from, unsigned class_to, double *damage, double *Sa_ret )
{
    double Sa = -1.0;  /* Negative amplitude indicates unset value */
    double D  =  0.0;

    assert( damage );
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut && !rfc_ctx->damage_lut_inapt)
    {
        if( !damage_calc_fast( rfc_ctx, class_from, class_to, &D, &Sa ) )
        {
            return false;
        }
    } else
#endif /*RFC_DAMAGE_FAST*/

#if RFC_USE_DELEGATES
    if( rfc_ctx->damage_calc_fcn )
    {
        if( !rfc_ctx->damage_calc_fcn( rfc_ctx, class_from, class_to, &D, &Sa ) )
        {
            return false;
        }
    } else
#endif /*RFC_USE_DELEGATES*/


    if( class_from != class_to )
    {
#if RFC_MINIMAL
        Sa = fabs( (int)class_from - (int)class_to ) / 2.0 * rfc_ctx->class_width;

        if( !damage_calc_amplitude( rfc_ctx, Sa, &D ) )
        {
            return false;
        }
#else /*!RFC_MINIMAL*/
        double Sa_i = fabs( (int)class_from - (int)class_to ) / 2.0 * rfc_ctx->class_width;
        double Sm_i =     ( (int)class_from + (int)class_to ) / 2.0 * rfc_ctx->class_width + rfc_ctx->class_offset;

        if( Sa_i > 0.0 )
        {
#if RFC_AT_SUPPORT
            /* Calculate transformation factor with normalized mean value */
            if( !RFC_at_transform( rfc_ctx, Sa_i, Sm_i, &Sa ) )
            {
                return false;
            }
#else /*!RFC_AT_SUPPORT*/
            Sa = Sa_i;
#endif /*RFC_AT_SUPPORT*/

            if( !damage_calc_amplitude( rfc_ctx, Sa, &D ) )
            {
                return false;
            }
        }
#endif /*RFC_MINIMAL*/
    }

    if( Sa_ret )
    {
        *Sa_ret = Sa;
    }

    *damage = D;

    return true;
}


#if RFC_DAMAGE_FAST
/**
 * @brief      Initialize a look-up table of damages for closed cycles. In this
 *             implementation the midrange doesn't matter!
 *
 * @param      rfc_ctx  The rainflow context
 */
static 
void damage_lut_init( rfc_ctx_s *rfc_ctx )
{
    double   *lut;
    unsigned  from, 
              to;
    double    Sa;

    assert( rfc_ctx );
    assert( rfc_ctx->damage_lut );
    assert( rfc_ctx->state == RFC_STATE_INIT );

    if( rfc_ctx->damage_lut )
    {
        lut = rfc_ctx->damage_lut;
        rfc_ctx->damage_lut = NULL;

        for( from = 0; from < rfc_ctx->class_count; from++ )
        {
            for( to = 0; to < rfc_ctx->class_count; to++ )
            {
                double D;

                if( !damage_calc( rfc_ctx, from, to, &D, &Sa ) )
                {
                    return;
                }
                lut[from * rfc_ctx->class_count + to] = D;
#if RFC_AT_SUPPORT
                if( rfc_ctx->amplitude_lut )
                {
                    rfc_ctx->amplitude_lut[from * rfc_ctx->class_count + to] = Sa;
                }
#endif /*RFC_AT_SUPPORT*/
            }
        }

        rfc_ctx->damage_lut       = lut;
        rfc_ctx->damage_lut_inapt = 0;
    }
}


/**
 * @brief      Calculate fictive damage for one closed (full) cycle, using
 *             look-up table.
 *
 * @param      rfc_ctx     The rainflow context
 * @param      class_from  The starting class
 * @param      class_to    The ending class
 * @param[out] damage      The damage value for the closed cycle
 * @param[out] Sa_ret      The amplitude (-1 if not available), may be NULL
 *
 * @return     true on success
 */
static
bool damage_calc_fast( rfc_ctx_s *rfc_ctx, unsigned class_from, unsigned class_to, double *damage, double *Sa_ret )
{
    double Sa = -1;    /* -1 means "uninitialized" */
    double D  =  0.0;

    assert( damage );
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );
    assert( class_from < rfc_ctx->class_count );
    assert( class_to   < rfc_ctx->class_count );

    if( rfc_ctx->damage_lut && !rfc_ctx->damage_lut_inapt )
    {
        D = rfc_ctx->damage_lut[class_from * rfc_ctx->class_count + class_to];

        if( Sa_ret )
        {
#if RFC_AT_SUPPORT
            if( rfc_ctx->amplitude_lut )
            {
                *Sa_ret = rfc_ctx->amplitude_lut[class_from * rfc_ctx->class_count + class_to];
            }
#else /*!RFC_AT_SUPPORT*/
            *Sa_ret = AMPLITUDE( rfc_ctx, fabs( (int)class_from - (int)class_to ) );
#endif /*RFC_AT_SUPPORT*/
        }
    }
    else
    {
        return error_raise( rfc_ctx, RFC_ERROR_LUT );
    }

    *damage = D;

    return true;
}
#endif /*RFC_DAMAGE_FAST*/


/**
 * @brief      Test data sample for a new turning point and add to the residue
 *             in that case. Update extrema.
 *             - 1. Hysteresis Filtering
 *             - 2. Peak-Valley Filtering
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  pt       The data point, must not be NULL
 *
 * @return     Returns pointer to new turning point in residue or NULL
 */
static
rfc_value_tuple_s * feed_filter_pt( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *pt )
{
    int                 slope;
    rfc_value_t         delta;
    rfc_value_tuple_s  *new_tp    = NULL;
    int                 do_append = 0;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state <= RFC_STATE_BUSY_INTERIM );

#if RFC_USE_DELEGATES
    if( rfc_ctx->tp_next_fcn )
    {
        return rfc_ctx->tp_next_fcn( rfc_ctx, pt );
    }
#endif /*RFC_USE_DELEGATES*/

    if( !pt ) return NULL;

    slope = rfc_ctx->internal.slope;

    /* Handle first turning point(s) */
    if( rfc_ctx->state < RFC_STATE_BUSY_INTERIM )
    {
        /* Residue is empty, still searching first turning point(s) */

        if( rfc_ctx->state == RFC_STATE_INIT )
        {
            /* Very first point, initialize local min-max search */
            rfc_ctx->internal.extrema[0] = 
            rfc_ctx->internal.extrema[1] = *pt;
            rfc_ctx->state               =  RFC_STATE_BUSY;
        }
        else
        {
            int is_falling_slope = -1;

            assert( rfc_ctx->state == RFC_STATE_BUSY );

            /* Still searching for first turning point(s) */

            /* Update local extrema */
            if( pt->value < rfc_ctx->internal.extrema[0].value )
            {
                /* Minimum */
                is_falling_slope = 1;
                rfc_ctx->internal.extrema[0] = *pt;
#if RFC_GLOBAL_EXTREMA
                rfc_ctx->internal.extrema_changed = true;
#endif /*RFC_GLOBAL_EXTREMA*/
            }
            else if( pt->value > rfc_ctx->internal.extrema[1].value )
            {
                /* Maximum */
                is_falling_slope = 0;
                rfc_ctx->internal.extrema[1] = *pt;
#if RFC_GLOBAL_EXTREMA
                rfc_ctx->internal.extrema_changed = true;
#endif /*RFC_GLOBAL_EXTREMA*/
            }

            /* Local hysteresis filtering */
            delta = value_delta( rfc_ctx->internal.extrema[0].value, rfc_ctx->internal.extrema[1].value, NULL /* sign_ptr */ );

            if( is_falling_slope >= 0 && delta > rfc_ctx->hysteresis )
            {
                /* Criteria met, new turning point found.
                 * Emit maximum on falling slope as first interim turning point, 
                 * minimum as second then (and vice versa) 
                 * 1st point: internal.extrema[ is_falling_slope]
                 * 2nd point: internal.extrema[!is_falling_slope]  ==> which is *pt also
                 */
                assert( rfc_ctx->residue_cnt < rfc_ctx->residue_cap );
                rfc_ctx->residue[rfc_ctx->residue_cnt] = rfc_ctx->internal.extrema[is_falling_slope];

                rfc_ctx->internal.slope = is_falling_slope ? -1 : 1;

                /* pt is the new interim turning point */
                rfc_ctx->state = RFC_STATE_BUSY_INTERIM;
                do_append = 1;
            }
        }
    }
    else  /* if( rfc_ctx->state < RFC_STATE_BUSY_INTERIM ) */
    {
        assert( rfc_ctx->state == RFC_STATE_BUSY_INTERIM );

        /* Consecutive search for turning points */

#if RFC_GLOBAL_EXTREMA
        /* Build global extrema */
        if( pt->value < rfc_ctx->internal.extrema[0].value )
        {
            /* Minimum */
            rfc_ctx->internal.extrema[0] = *pt;
            rfc_ctx->internal.extrema_changed = true;
        }
        else if( pt->value > rfc_ctx->internal.extrema[1].value )
        {
            /* Maximum */
            rfc_ctx->internal.extrema[1] = *pt;
            rfc_ctx->internal.extrema_changed = true;
        }
#endif /*RFC_GLOBAL_EXTREMA*/

        /* Hysteresis Filtering, check against interim turning point */
        delta = value_delta( rfc_ctx->residue[rfc_ctx->residue_cnt].value, pt->value, &slope /* sign_ptr */ );

        /* There are three scenarios possible here:
         *   1. Previous slope is continued
         *      "delta" is ignored whilst hysteresis is exceeded already.
         *      Interim turning point has just to be adjusted.
         *   2. Slope reversal, slope is greater than hysteresis
         *      Interim turning point becomes real turning point.
         *      Current point becomes new interim turning point
         *   3. Slope reversal, slope is less than or equal hysteresis
         *      Nothing to do.
         */

        /* Peak-Valley Filtering */
        /* Justify interim turning point, or add a new one (2) */
        if( slope == rfc_ctx->internal.slope )
        {
            /* Scenario (1), Continuous slope */

            /* Replace interim turning point with new extrema */
            if( rfc_ctx->residue[rfc_ctx->residue_cnt].value != pt->value )
            {
                rfc_ctx->residue[rfc_ctx->residue_cnt] = *pt;
            }
        }
        else
        {
            if( delta > rfc_ctx->hysteresis )
            {
                /* Scenario (2), Criteria met: slope != rfc_ctx->internal.slope && delta > rfc_ctx->hysteresis */

                /* Storage */
                rfc_ctx->internal.slope = slope;

                /* Handle new turning point */
                do_append = 1;
            }
            else
            {
                /* Scenario (3), Turning point found, but still in hysteresis band, nothing to do */
            }
        }
    }

    /* Handle new turning point, that is the current last point in residue */
    if( do_append )
    {
        assert( rfc_ctx->state == RFC_STATE_BUSY_INTERIM );

        /* Increment and set new interim turning point */
        assert( rfc_ctx->residue_cnt < rfc_ctx->residue_cap );
        rfc_ctx->residue[++rfc_ctx->residue_cnt] = *pt;

        /* Return new turning point */
        new_tp = &rfc_ctx->residue[rfc_ctx->residue_cnt - 1];
    }

    return new_tp;
}


# if !RFC_MINIMAL
/**
 * @brief      Rainflow counting core, assumes
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 */
static
void cycle_find( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->cycle_find_fcn && rfc_ctx->counting_method == RFC_COUNTING_METHOD_DELEGATED )
    {
        rfc_ctx->cycle_find_fcn( rfc_ctx, flags );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        switch( rfc_ctx->counting_method )
        {
            case RFC_COUNTING_METHOD_NONE:
                break;
            case RFC_COUNTING_METHOD_4PTM:
                cycle_find_4ptm( rfc_ctx, flags );
                break;
#if RFC_HCM_SUPPORT
            case RFC_COUNTING_METHOD_HCM:
                cycle_find_hcm( rfc_ctx, flags );
                break;
#endif /*RFC_HCM_SUPPORT*/
            case RFC_COUNTING_METHOD_DELEGATED:
                /* FALLTHROUGH */
            default:
                assert( false );
                break;
        }
    }

    /* If no rainflow counting is done, just look out for turning points, discard residue */
    if( rfc_ctx->counting_method == RFC_COUNTING_METHOD_NONE || !rfc_ctx->class_count )
    {
        /* Prune residue */
        if( rfc_ctx->residue_cnt > 1 )
        {
            residue_remove_item( rfc_ctx, /*index*/ 0, rfc_ctx->residue_cnt - 1 );
        }
    }

}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Rainflow counting core (4-point-method [3]).
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 */
static
void cycle_find_4ptm( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    while( rfc_ctx->residue_cnt >= 4 )
    {
        size_t idx = rfc_ctx->residue_cnt - 4;

        unsigned A = rfc_ctx->residue[idx+0].cls;
        unsigned B = rfc_ctx->residue[idx+1].cls;
        unsigned C = rfc_ctx->residue[idx+2].cls;
        unsigned D = rfc_ctx->residue[idx+3].cls;

        if( B > C )
        {
            unsigned temp = B;
            B = C;
            C = temp;
        }

        if( A > D )
        {
            unsigned temp = A;
            A = D;
            D = temp;
        }

        /* Check for closed cycles [3] */
        if( A <= B && C <= D )
        {
            rfc_value_tuple_s *from = &rfc_ctx->residue[idx+1];
            rfc_value_tuple_s *to   = &rfc_ctx->residue[idx+2];

            /* Closed cycle found, process countings */
            cycle_process_counts( rfc_ctx, from, to, to + 1, flags );

            /* Remove two inner turning points (idx+1 and idx+2) */
            /* Move last turning point */
            rfc_ctx->residue[idx+1] = rfc_ctx->residue[idx+3];
            /* Move interim turning point */
            if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
            {
                rfc_ctx->residue[idx+2] = rfc_ctx->residue[idx+4];
            }
            rfc_ctx->residue_cnt -= 2;
        }
        else break;
    }
}


#if RFC_HCM_SUPPORT
/**
 * @brief      Rainflow counting core (HCM method).
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 */
static
void cycle_find_hcm( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    int IZ, IR;
    double eps = rfc_ctx->class_width / 100;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    IZ = rfc_ctx->internal.hcm.IZ - 1,  /* hcm.IZ and hcm.IR are base 1! */
    IR = rfc_ctx->internal.hcm.IR - 1;

    while( rfc_ctx->residue_cnt > 0 )
    {
        rfc_value_tuple_s *I, *J, *K;

        /* Translation from "RAINFLOW.F" */
/*label_1:*/
        K = rfc_ctx->residue;  /* Recent value (turning point) */

        /* Place first turning point into stack */
        if( !IR )
        {
            rfc_ctx->internal.hcm.stack[IR++] = *K;
        }

label_2:
        if( IZ > IR )
        {
            /* There are at least 2 cycles on the stack able to close */
            I = &rfc_ctx->internal.hcm.stack[IZ-1];
            J = &rfc_ctx->internal.hcm.stack[IZ];

            if( (double)(K->value - J->value) * (double)(J->value - I->value) + eps >= 0.0 )
            {
                /* Is no turning point */
                /* This should only may happen, when RFC_FLAGS_ENFORCE_MARGIN is set, 
                   since all values from residue are turning points */
                assert( rfc_ctx->internal.flags & RFC_FLAGS_ENFORCE_MARGIN );
                IZ--;
                /* Test further closed cycles */
                goto label_2;
            }
            else
            {
                /* Is a turning point */
                if( fabs( (double)K->value - (double)J->value ) + eps >= fabs( (double)J->value - (double)I->value ) )
                {
                    /* Cycle range is greater or equal to previous, register closed cycle */
                    cycle_process_counts( rfc_ctx, I, J, NULL, flags );
                    IZ -= 2;
                    /* Test further closed cycles */
                    goto label_2;
                }
            }
        }
        else if( IZ == IR )
        {
            J = &rfc_ctx->internal.hcm.stack[IZ];

            if( ( (double)K->value - (double)J->value ) * (double)J->value + eps >= 0.0 )
            {
                /* Is no turning point */
                IZ--;
                /* Test further closed cycles */
                goto label_2;
            }
            else if( fabs( (double)K->value ) + eps > fabs( (double)J->value ) )
            {
                /* Is turning point and range is less than previous */
                IR++;
            }
        }
        else
        {
            /* IZ < IR: There is no cycle on the stack able to close */
        }

        /* Place cycle able to close */
        IZ++;
        assert( IZ < rfc_ctx->internal.hcm.stack_cap );
        rfc_ctx->internal.hcm.stack[IZ] = *K;

        /* "goto" not necessary: while loop */
        /* goto label_1; */

        /* Remove K (first element) from residue */
        residue_remove_item( rfc_ctx, /*index*/ 0, 1 );
    }

    /* hcm.IZ and hcm.IR are base 1! */
    rfc_ctx->internal.hcm.IZ = IZ + 1;
    rfc_ctx->internal.hcm.IR = IR + 1;
}
#endif /*RFC_HCM_SUPPORT*/


#if !RFC_MINIMAL
/**
 * @brief      Processes LC count (level crossing) for the last slope in residue
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    Control flags
 */
static
void cycle_process_lc( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    size_t n = rfc_ctx->residue_cnt;

    if( n > 1 && (flags & RFC_FLAGS_COUNT_LC) )
    {
        /* Do the level crossing counting */
        bool rising = rfc_ctx->residue[n-1].value > rfc_ctx->residue[n-2].value;

        if( rising )
        {
            cycle_process_counts( rfc_ctx, &rfc_ctx->residue[n-2], &rfc_ctx->residue[n-1], NULL, flags & (RFC_FLAGS_COUNT_LC_UP | RFC_FLAGS_ENFORCE_MARGIN) );
        }
        else
        {
            cycle_process_counts( rfc_ctx, &rfc_ctx->residue[n-2], &rfc_ctx->residue[n-1], NULL, flags & (RFC_FLAGS_COUNT_LC_DN | RFC_FLAGS_ENFORCE_MARGIN) );
        }
    }
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief         Processes counts on a closing cycle and further modifications
 *                due to the new damage value.
 *
 * @param         rfc_ctx  The rainflow context
 * @param[in,out] from     The starting data point
 * @param[in,out] to       The ending data point
 * @param[in,out] next     The point next after "to"
 * @param         flags    Control flags
 */
static
void cycle_process_counts( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, rfc_flags_e flags )
{
    unsigned class_from, class_to;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );
    assert( !rfc_ctx->class_count || ( from->value > rfc_ctx->class_offset && to->value > rfc_ctx->class_offset ) );

#if RFC_TP_SUPPORT
    /* If flag RFC_FLAGS_ENFORCE_MARGIN is set, cycles less than hysteresis are possible */
    if( flags & RFC_FLAGS_ENFORCE_MARGIN )
    {
        if( value_delta( from->value, to->value, NULL /* sign_ptr */ ) <= rfc_ctx->hysteresis )
        {
            return;
        }
    }
#endif /*RFC_TP_SUPPORT*/

    /* Quantized "from" */
    class_from = from->cls;

    if( class_from >= rfc_ctx->class_count ) class_from = rfc_ctx->class_count - 1;

    /* Quantized "to" */
    class_to = to->cls;

    if( class_to >= rfc_ctx->class_count ) class_to = rfc_ctx->class_count - 1;
    
    /* class_from and class_to are base 0 now */

    /* Do several counts, according to "flags" */
    if( class_from != class_to )
    {
#if RFC_DEBUG_FLAGS
        if( rfc_ctx->internal.debug_flags & RFC_FLAGS_LOG_CLOSED_CYCLES &&
            flags & (RFC_FLAGS_COUNT_ALL & ~RFC_FLAGS_COUNT_LC) )
        {
            RFC_debug_fprintf( rfc_ctx, stdout, 
                               "Closed cycle %g @ %lu[%lu] --> %g @ %lu[%lu]\n", 
                               from->value, (long unsigned)from->pos, (long unsigned)from->tp_pos,
                               to->value,   (long unsigned)to->pos,   (long unsigned)to->tp_pos );
        }
#endif /*RFC_DEBUG_FLAGS*/
        /* Cumulate damage */
        if( flags & RFC_FLAGS_COUNT_DAMAGE )
        {
            double Sa_i;
            double D_i;

            if( !damage_calc( rfc_ctx, class_from, class_to, &D_i, &Sa_i ) )
            {
                return;
            }

            /* Adding damage for the current cycle, with its actual weight */
            rfc_ctx->damage += D_i * rfc_ctx->curr_inc / rfc_ctx->full_inc;
#if !RFC_MINIMAL
            /* Fatigue strength Sd(D) depresses in subject to cumulative damage D.
               Sd(D)/Sd = (1-D)^(1/q), [6] chapter 3.2.9, formula 3.2-44 and 3.2-46
               Only cycles exceeding Sd(D) have damaging effect. */
            if( Sa_i >= rfc_ctx->internal.wl.sd && ( flags & RFC_FLAGS_COUNT_MK ) )
            {
                rfc_wl_param_s  wl_unimp;                          /* WL parameters unimpaired part */
                rfc_wl_param_s *wl_imp = &rfc_ctx->internal.wl;    /* WL parameters impaired part */
                double          D_con;                             /* Current damage, Miners' consequent rule */

                /* Backup Woehler curve parameters and use shadowed ones for the impaired part instead */
                RFC_wl_param_get( rfc_ctx, &wl_unimp );
                RFC_wl_param_set( rfc_ctx,  wl_imp );

#if RFC_DAMAGE_FAST
                if( rfc_ctx->damage_lut )
                {
                    /* Disable lut temporarily, since it is only valid for Woehler parameters for unimpaired part */
                    rfc_ctx->damage_lut_inapt++;
                    (void)damage_calc_amplitude( rfc_ctx, Sa_i, &D_con );
                    rfc_ctx->damage_lut_inapt--;
                }
                else
#endif /*RFC_DAMAGE_FAST*/
                {
                    (void)damage_calc_amplitude( rfc_ctx, Sa_i, &D_con );
                }

                D_con += wl_imp->D;

                if( D_con < 1.0 )
                {
                    /* Calculate new parameters for the Woehler curve, impaired part */
                    if( wl_unimp.sx > 0.0 )
                    {
                        double q       =       wl_imp->q;
                        double k       = fabs( wl_imp->k );
                        rfc_ctx->wl_sx = wl_unimp.sx * pow( 1.0 - D_con, 1.0 / q );
                     /* rfc_ctx->wl_nx = wl.nx * pow( 1.0 - D_con, -( k - q ) / q );  // Damage accumulation is done in D_con! */

                        (void)RFC_wl_calc_n( rfc_ctx, 
                                             wl_unimp.sx, 
                                             wl_unimp.nx, 
                                             k,
                                             rfc_ctx->wl_sx,
                                            &rfc_ctx->wl_nx );
                    }

                    if( wl_unimp.sd > 0.0 )
                    {
                        double q2      =       wl_imp->q2;
                        double k2      = fabs( wl_imp->k2 );
                        rfc_ctx->wl_sd = wl_unimp.sd * pow( 1.0 - D_con, 1.0 / q2 );
                     /* rfc_ctx->wl_nx = wl_unimp.nx * pow( 1.0 - D_con, -( k2 - q2 ) / q2 );  // Damage accumulation is done in D_con! */

                        (void)RFC_wl_calc_n( rfc_ctx, 
                                             wl_unimp.sd, 
                                             wl_unimp.nd, 
                                             k2,
                                             rfc_ctx->wl_sd,
                                            &rfc_ctx->wl_nd );
                    }
                }

                RFC_wl_param_get( rfc_ctx, wl_imp );
                rfc_ctx->internal.wl.D = D_con;
                RFC_wl_param_set( rfc_ctx, &wl_unimp );
            }
#endif /*!RFC_MINIMAL*/
        }

        /* Rainflow matrix */
        if( rfc_ctx->rfm && ( flags & RFC_FLAGS_COUNT_RFM ) )
        {
            /* Matrix (row-major storage):
             *          t o
             *    +-------------
             *    | 0 1 2 3 4 5
             *  f | 6 7 8 9 # #
             *  r | # # # # # #
             *  o | # # # # # #
             *  m | # # # # # #
             *    | # # # # # #<-(6x6-1)
             */
            size_t idx = rfc_ctx->class_count * class_from + class_to;
            
            assert( rfc_ctx->rfm[idx] <= RFC_COUNTS_LIMIT );
            rfc_ctx->rfm[idx] += rfc_ctx->curr_inc;
        }

#if !RFC_MINIMAL
        /* Range pair */
        if( rfc_ctx->rp && ( flags & RFC_FLAGS_COUNT_RP ) )
        {
            /* 
             * Range pair histogram (vector storage)
             * Range value = idx * class_width  (=2x Amplitude)
             */
            int idx = abs( (int)class_from - (int)class_to );
            
            assert( rfc_ctx->rp[idx] <= RFC_COUNTS_LIMIT );
            rfc_ctx->rp[idx] += rfc_ctx->curr_inc;
        }

        /* Level crossing, count rising and falling slopes */
        if( rfc_ctx->lc && ( flags & RFC_FLAGS_COUNT_LC ) )
        {
            /* 
             * Level crossing histogram (vector storage)
             * Counts class upper bound crossings
             * Class upper bound value = (idx+1) * class_width + class_offset
             */
            unsigned idx;
            unsigned idx_from = ( class_from < class_to ) ? class_from : class_to;
            unsigned idx_to   = ( class_from > class_to ) ? class_from : class_to;

            for( idx = idx_from; idx < idx_to; idx++ )
            {
                if( flags & RFC_FLAGS_COUNT_LC_UP )
                {
                    /* Count rising slopes */
                    assert( rfc_ctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfc_ctx->lc[idx] += rfc_ctx->full_inc;
                }

                if( flags & RFC_FLAGS_COUNT_LC_DN )
                {
                    /* Count falling slopes */
                    assert( rfc_ctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfc_ctx->lc[idx] += rfc_ctx->full_inc;
                }
            }
        }

#if RFC_DH_SUPPORT
        if( flags & RFC_FLAGS_COUNT_DH )
        {
            /* "Spread" damage over turning points (tp) and damage history (dh) */
            spread_damage( rfc_ctx, from, to, next, flags );
        }
#endif /*RFC_DH_SUPPORT*/

#endif /*!RFC_MINIMAL*/
    }
}


#if RFC_TP_SUPPORT
/**
 * @brief         Append or alter a turning point in its storage.
 *
 *                Attention: Consider tp_locked! tp_cnt and tp_cap have to
 *                reflect storage state!
 *
 * @param         rfc_ctx    The rainflow context
 * @param[in]     tp_pos     The position. If tp_pos==0 and tp->tp_pos==0, tp
 *                           will be appended. If tp_pos==0 and tp->tp_pos>0,
 *                           only the position tp->tp_pos is altered. If
 *                           tp_pos>0 the existing turning point is overwritten
 *                           by tp.
 * @param[in,out] tp         The turning point
 *
 * @return        true on success
 */
static
bool tp_set( rfc_ctx_s *rfc_ctx, size_t tp_pos, rfc_value_tuple_s *tp )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state <= RFC_STATE_FINISHED );
    assert( tp );

    if( !tp || rfc_ctx->tp_locked )
    {
        return false;
    }

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->tp_set_fcn )
    {
        /* Add turning point */
        return rfc_ctx->tp_set_fcn( rfc_ctx, tp_pos, tp );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        if( !rfc_ctx->tp )
        {
            return true;
        }

        /* Check to append or alter */
        if( tp_pos )
        {
            /* Alter or move existing turning point */
            assert( tp_pos <= rfc_ctx->tp_cnt );

            tp->tp_pos                =  0;        /* No position information for turning points in its storage */
            rfc_ctx->tp[ tp_pos - 1 ] = *tp;       /* Move or replace turning point */
            tp->tp_pos                =  tp_pos;   /* Ping back the position (commonly tp lies in residue buffer) */

            return true;
        }
        else
        {
            /* Append (tp_pos == 0) */
            if( tp->tp_pos )
            {
                /* Already an element of tp stack */
                assert( tp->tp_pos <= rfc_ctx->tp_cnt );

                return true;
            }
            else
            {
                /* Append tp at the tail */
                tp_pos = ++rfc_ctx->tp_cnt;
            }
        }

        /* Check if buffer needs to be resized */
        if( rfc_ctx->tp_cnt >= rfc_ctx->tp_cap )
        {
            rfc_value_tuple_s  *tp_new;
            size_t              tp_cap_new;
            size_t              tp_cap_increment;

            /* Reallocation */
            tp_cap_increment = (size_t)1024 * ( rfc_ctx->tp_cap / 640 + 1 );  /* + 60% + 1024 */
            tp_cap_new       = rfc_ctx->tp_cap + tp_cap_increment;
            tp_new           = rfc_ctx->mem_alloc( rfc_ctx->tp, tp_cap_new, 
                                                   sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TP );

            if( tp_new )
            {
                rfc_ctx->tp     = tp_new;
                rfc_ctx->tp_cap = tp_cap_new;
            }
            else
            {
                return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
            }
        }

        assert( tp_pos <= rfc_ctx->tp_cnt );

        /* Append turning point */
        rfc_ctx->tp[ tp_pos - 1 ] = *tp;      /* Make a copy of tp in .tp, tp->tp_pos remains unaltered */
        tp->tp_pos                =  tp_pos;  /* Ping back turning point position index in tp, base 1 */

        if( rfc_ctx->internal.flags & RFC_FLAGS_TPAUTOPRUNE && rfc_ctx->tp_cnt > rfc_ctx->tp_prune_threshold )
        {
            return RFC_tp_prune( rfc_ctx, rfc_ctx->tp_prune_size, RFC_FLAGS_TPPRUNE_PRESERVE_POS );
        }

        return true;
    }
}


/**
 * @brief      Get turning point reference
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  tp_pos   The position, base 1
 * @param      tp       The buffer for the reference
 *
 * @return     true on success
 */
static
bool tp_get( rfc_ctx_s *rfc_ctx, size_t tp_pos, rfc_value_tuple_s **tp )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state <= RFC_STATE_FINISHED );

    if( !tp_pos || tp_pos > rfc_ctx->tp_cnt )
    {
        return false;
    }

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->tp_get_fcn )
    {
        /* Add turning point */
        return rfc_ctx->tp_get_fcn( rfc_ctx, tp_pos, tp );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        if( !rfc_ctx->tp )
        {
            return false;
        }

        if( tp )
        {
            assert( tp_pos <= rfc_ctx->tp_cnt );
            *tp = &rfc_ctx->tp[ tp_pos - 1 ];
        }
    }

    return true;
}


/**
 * @brief      Increase damage for existing turning point.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  tp_pos   The tp position (base 1)
 * @param[in]  damage   The damage
 *
 * @return     true on success
 */
static
bool tp_inc_damage( rfc_ctx_s *rfc_ctx, size_t tp_pos, double damage )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    /* Add new turning point */

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->tp_inc_damage_fcn )
    {
        /* Add turning point */
        return rfc_ctx->tp_inc_damage_fcn( rfc_ctx, tp_pos, damage );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        if( rfc_ctx->tp && tp_pos )
        {
            if( !tp_pos || tp_pos > rfc_ctx->tp_cnt )
            {
                return error_raise( rfc_ctx, RFC_ERROR_TP );
            }
#if RFC_DH_SUPPORT
            rfc_ctx->tp[ tp_pos - 1 ].damage += damage;
#endif /*RFC_DH_SUPPORT*/
        }
    }

    return true;
}


/**
 * @brief      Lock turning points queue
 *
 * @param      rfc_ctx  The rainflow context
 * @param      do_lock  Turning point storage will be locked, if true
 */
static 
void tp_lock( rfc_ctx_s *rfc_ctx, bool do_lock )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    rfc_ctx->tp_locked += do_lock ? 1 : -1;

    if( rfc_ctx->tp_locked < 0 )
    {
        rfc_ctx->tp_locked = 0;
    }
}


/**
 * @brief      Restart counting with given points from turning points history
 *
 * @param      rfc_ctx          The rainflow context
 * @param      new_hysteresis   The new hysteresis
 * @param[in]  new_class_param  The new class parameters
 *
 * @return     true on success
 */
static
bool tp_refeed( rfc_ctx_s *rfc_ctx, rfc_value_t new_hysteresis, const rfc_class_param_s *new_class_param )
{
    size_t pos_offset,
           i;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    if( !feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    pos_offset = rfc_ctx->internal.pos_offset;

    /* Clear current count data */
    RFC_clear_counts( rfc_ctx );

    rfc_ctx->internal.pos_offset = pos_offset;

    /* Class parameters may change, new hysteresis must be greater! */
    if( new_class_param )
    {
        assert( new_hysteresis >= rfc_ctx->hysteresis );
        rfc_ctx->hysteresis     = new_hysteresis;
#if RFC_DAMAGE_FAST
        rfc_ctx->damage_lut_inapt = 1;
        if( rfc_ctx->class_count != new_class_param->count )
        {
            size_t num = rfc_ctx->class_count * rfc_ctx->class_count;
            if( !rfc_ctx->mem_alloc( rfc_ctx->damage_lut, num, sizeof(rfc_counts_t), RFC_MEM_AIM_DLUT ) )
            {
                return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
            }
#if RFC_AT_SUPPORT
            if( rfc_ctx->amplitude_lut )
            {
                if( !rfc_ctx->mem_alloc( rfc_ctx->amplitude_lut, num, sizeof(rfc_counts_t), RFC_MEM_AIM_DLUT ) )
                {
                    return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
                }
            }
#endif /*RFC_AT_SUPPORT*/
        }
        RFC_class_param_set( rfc_ctx, new_class_param );
        damage_lut_init( rfc_ctx );
#else /*!RFC_DAMAGE_FAST*/
        RFC_class_param_set( rfc_ctx, new_class_param );
#endif /*RFC_DAMAGE_FAST*/
    }

#if RFC_USE_DELEGATES
    if( rfc_ctx->tp_get_fcn )
    {
        rfc_value_tuple_s *tp, *tp_ext;
        size_t             tp_cnt = rfc_ctx->tp_cnt;
        int                n      = 0;
        const int          n_max  = 500;
        bool               ok     = true;  

        tp = mem_alloc( NULL, n_max, sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TP );

        if( !tp )
        {
            return error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }

        rfc_ctx->tp_cnt = 0;
        i = 1;
        while( tp_cnt && ok )
        {
            int n_cnt = ( tp_cnt <= (size_t)n_max ) ? (int)tp_cnt : n_max;

            for( n = 0; ok && n < n_cnt; i++, n++ )
            {
                if( !tp_get( rfc_ctx, i, &tp_ext ) )
                {
                    ok = false;
                    break;
                }
                else
                {
                    tp[n]        = *tp_ext;
                    tp[n].cls    =  QUANTIZE( rfc_ctx, tp[n].value );
#if RFC_DH_SUPPORT
                    tp[n].damage = 0.0;
#endif /*RFC_DH_SUPPORT*/
                }
            }

            RFC_feed_tuple( rfc_ctx, tp, n_cnt );

            tp_cnt -= n_cnt;
        }

        mem_alloc( tp, 0, 0, RFC_MEM_AIM_TP );

        return ok;
    }
    else
#endif /*!RFC_USE_DELEGATES*/
    {
        size_t             tp_cnt = rfc_ctx->tp_cnt;
        rfc_value_tuple_s *tp     = rfc_ctx->tp;

        for( i = 0; i < tp_cnt; i++, tp++ )
        {
            tp->cls    = QUANTIZE( rfc_ctx, tp->value );
#if RFC_DH_SUPPORT
            tp->damage = 0.0;
#endif /*RFC_DH_SUPPORT*/
        }

        rfc_ctx->tp_cnt = 0;
        RFC_feed_tuple( rfc_ctx, tp, tp_cnt );
    }
}
#endif /*RFC_TP_SUPPORT*/


#if RFC_DH_SUPPORT
/**
 * @brief         Spread damage over turning points and damage history
 *
 * @param         rfc_ctx  The rainflow context
 * @param[in,out] from     The start turning point
 * @param[in,out] to       The end turning point
 * @param[in,out] next     The next turning point after @a to
 * @param         flags    The flags
 *
 * @return        true on success
 */
static 
bool spread_damage( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *from, 
                                        rfc_value_tuple_s *to, 
                                        rfc_value_tuple_s *next, rfc_flags_e flags )
{
    int    spread_damage_method;
    double D = 0.0;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );
    assert( from && to );

    spread_damage_method = rfc_ctx->spread_damage_method;

#if RFC_TP_SUPPORT
    if( !from->tp_pos && !to->tp_pos )
    {
        return true;
    }

    if( !from->tp_pos ) spread_damage_method = RFC_SD_FULL_P3;
    if( !to->tp_pos )   spread_damage_method = RFC_SD_FULL_P2;
#endif /*RFC_TP_SUPPORT*/

    switch( spread_damage_method  )
    {
        case RFC_SD_NONE:
            break;
        case RFC_SD_HALF_23:
        case RFC_SD_FULL_P2:
        case RFC_SD_FULL_P3:
        {
            double damage_lhs, damage_rhs;

            if( !damage_calc( rfc_ctx, from->cls, to->cls, &D, NULL /*Sa_ret*/ ) )
            {
                return false;
            }

            /* Current cycle weight */
            D *= (double)rfc_ctx->curr_inc / rfc_ctx->full_inc;

            if( rfc_ctx->spread_damage_method == RFC_SD_FULL_P2 )
            {
                damage_lhs = D;
                damage_rhs = 0.0;
            }
            else if( rfc_ctx->spread_damage_method == RFC_SD_FULL_P3 )
            {
                damage_lhs = 0.0;
                damage_rhs = D;
            }
            else
            {
                damage_lhs = damage_rhs = D / 2.0;
            }

#if RFC_TP_SUPPORT
            if( from->tp_pos && !tp_inc_damage( rfc_ctx, from->tp_pos, damage_lhs ) )
            {
                return false;
            }

            if( to->tp_pos && !tp_inc_damage( rfc_ctx, to->tp_pos, damage_rhs ) )
            {
                return false;
            }
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
            if( rfc_ctx->dh )
            {
                if( from->pos )
                {
                    rfc_ctx->dh[ from->pos - 1 ] += damage_lhs;
                }
                else damage_rhs += damage_lhs;

                if( to->pos )
                {
                    rfc_ctx->dh[ to->pos - 1 ] += damage_rhs;
                }
            }
#endif /*RFC_DH_SUPPORT*/

            break;
        }
        case RFC_SD_RAMP_AMPLITUDE_23:
        case RFC_SD_RAMP_DAMAGE_23:
        case RFC_SD_RAMP_AMPLITUDE_24:
        case RFC_SD_RAMP_DAMAGE_24:
        {
#if RFC_TP_SUPPORT
            size_t  i,
                    start, end,         /* Base 0 */
                    width, 
                    tp_start, tp_end;   /* Base 0 */
            int     from_cls, to_cls;   /* Base 0 */

            /* Care about possible wrapping, caused by RFC_RES_REPEATED:
                0         1
                01234567890123 (history of 14 turning points)
                ...E....S.....
                   ^End ^Start
                   =3   =8

                Results in:
                0         1         2
                0123456789012345678901234567
                ...E....S....,...E....S.....
                        ^Start   ^End
                        =8       =17
            */

            from_cls = from->cls;
            to_cls   = to->cls;

            /* Spread over P2 to P3 or over P2 to P4 */
            if( rfc_ctx->spread_damage_method == RFC_SD_RAMP_AMPLITUDE_24 ||
                rfc_ctx->spread_damage_method == RFC_SD_RAMP_DAMAGE_24 )
            {
                to = next ? next : to;
            }
            
            /* Absolute position (input stream) */
            start    = from->pos - 1;
            end      = to->pos - 1;
            end     += ( start >= end ) ? rfc_ctx->internal.pos : 0;   /* internal.pos is not modified while repeated counting */
            width    = end - start;
            /* Position in turning point storage */
            tp_start = from->tp_pos - 1;
            tp_end   = to->tp_pos - 1;
            tp_end  += ( tp_start >= tp_end ) ? rfc_ctx->tp_cnt : 0;   /* tp is modified while repeated counting, but tp_pos is consistent */

            D        = 0.0;

            assert( width > 0 );

            /* Iterate over turning points */
            for( i = tp_start; i <= tp_end; i++ )
            {
                rfc_value_tuple_s *tp;
                size_t             tp_pos_0,        /* Turning point position, base 0 */
                                   pos_0;           /* Input stream position, base 0 */
                double             weight;
                double             D, D_new;

                tp_pos_0 = i % rfc_ctx->tp_cnt;

                if( !tp_get( rfc_ctx, tp_pos_0 + 1, &tp ) )
                {
                    return false;
                }

                pos_0    = tp->pos - 1;
                pos_0   += ( start >= pos_0 ) ? rfc_ctx->internal.pos : 0;
                weight   = (double)( pos_0 - start ) / width;

                switch( rfc_ctx->spread_damage_method )
                {
                    case RFC_SD_RAMP_AMPLITUDE_23:
                    case RFC_SD_RAMP_AMPLITUDE_24:
                        if( !damage_calc( rfc_ctx, from_cls, to_cls, &D, NULL /*Sa_ret*/ ) )
                        {
                            return false;
                        }
                        /* Current cycle weight */
                        D *= (double)rfc_ctx->curr_inc / rfc_ctx->full_inc;
                        /* Di = 1/( ND*(Sa/SD*weight)^k ) = 1/( ND*(Sa/SD)^k ) * 1/weight^k */
                        D_new = D * pow( weight, -fabs(rfc_ctx->wl_k) );
                        break;
                    case RFC_SD_RAMP_DAMAGE_23:
                    case RFC_SD_RAMP_DAMAGE_24:
                        if( !damage_calc( rfc_ctx, from_cls, to_cls, &D, NULL /*Sa_ret*/ ) )
                        {
                            return false;
                        }
                        /* Current cycle weight */
                        D *= (double)rfc_ctx->curr_inc / rfc_ctx->full_inc;
                        /* Di = 1/( ND*(Sa/SD)^k ) * weight */
                        D_new = D * weight;
                        break;
                }

                if( D_new > D )
                {
#if RFC_DH_SUPPORT
                    if( rfc_ctx->dh )
                    {
                        rfc_ctx->dh[ pos_0 ] += D_new - D;
                    }
#endif /*RFC_DH_SUPPORT*/
                    if( !tp_inc_damage( rfc_ctx, tp_pos_0 + 1, D_new - D ) )
                    {
                        return false;
                    }
                    D = D_new;
                }
            }
#endif /*RFC_TP_SUPPORT*/
        }
        break;

        case RFC_SD_TRANSIENT_23:
            /* \todo */
            return error_raise( rfc_ctx, RFC_ERROR_UNSUPPORTED );
            break;
        case RFC_SD_TRANSIENT_23c:
            /* \todo */
            return error_raise( rfc_ctx, RFC_ERROR_UNSUPPORTED );
            break;

        default:
            assert( false );
            break;
    }

    return true;
}
#endif /*RFC_DH_SUPPORT*/


/**
 * @brief      Raises an error
 *
 * @param      rfc_ctx  The rainflow context
 * @param      error    The error identifier
 *
 * @return     false on error
 */
static
bool error_raise( rfc_ctx_s *rfc_ctx, rfc_error_e error )
{
    if( error == RFC_ERROR_NOERROR ) return true;

    if( rfc_ctx && rfc_ctx->version == sizeof(rfc_ctx_s) )
    {
        rfc_ctx->state = RFC_STATE_ERROR;
        rfc_ctx->error = error;
    }
    else
    {
        assert( false );
    }

    return false;
}


/**
 * @brief      Returns the unsigned difference of two values, sign optionally
 *             returned as -1 or 1.
 *
 * @param      from_val  Left hand value
 * @param      to_val    Right hand value
 * @param[out] sign_ptr  Pointer to catch sign (may be NULL)
 *
 * @return     Returns the absolute difference of given values
 */
static
rfc_value_t value_delta( rfc_value_t from_val, rfc_value_t to_val, int *sign_ptr )
{
    double delta = (double)to_val - (double)from_val;

    if( sign_ptr )
    {
        *sign_ptr = ( delta < 0.0 ) ? -1 : 1;
    }

    return (rfc_value_t)fabs( delta );
}


/**
 * @brief      (Re-)Allocate or free memory
 *
 * @param      ptr   Previous data pointer, or NULL, if unset
 * @param      num   The number of elements
 * @param      size  The size of one element in bytes
 * @param      aim   The aim
 *
 * @return     New memory pointer or NULL if either num or size is 0
 */
static
void * mem_alloc( void *ptr, size_t num, size_t size, rfc_mem_aim_e aim )
{
    if( !num || !size )
    {
        if( ptr )
        {
            free( ptr );
        }
        return NULL;
    }
    else
    {
        return ptr ? realloc( ptr, num * size ) : calloc( num, size );
    }
}




/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


#if MATLAB_MEX_FILE

#if RFC_DEBUG_FLAGS
static
int rfc_vfprintf_fcn( void *ctx, FILE* stream, const char *fmt, va_list arg )
{
    int length;

    /* Get the buffer size needed */
    length = vsnprintf( NULL, 0, fmt, arg );  
    if( length > 0 )
    {
        char *buffer = malloc( ++length );

        if( buffer )
        {
            buffer[length-1] = 0;
            vsnprintf( buffer, length, fmt, arg );
            mexPrintf( "%s", buffer );

            free( buffer );
        }
    }

    return length;
}
#endif /*RFC_DEBUG_FLAGS*/

/**
 * MATLAB wrapper for the rainflow algorithm
 */
static
void mexRainflow( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
#if !RFC_MINIMAL
    if( nrhs != 9 )
    {
        mexErrMsgTxt( "Function needs exact 8 arguments!" );
#else /*RFC_MINIMAL*/
    if( nrhs != 5 )
    {
        if( !nrhs )
        {
            mexPrintf( "%s", RFC_MEX_USAGE );
            return;
        }
        mexErrMsgTxt( "Function needs exact 5 arguments!" );

#endif /*!RFC_MINIMAL*/
    }
    else
    {
        rfc_ctx_s rfc_ctx = { sizeof(rfc_ctx_s) };
    
        const mxArray  *mxData           = prhs[0];
        const mxArray  *mxClassCount     = prhs[1];
        const mxArray  *mxClassWidth     = prhs[2];
        const mxArray  *mxClassOffset    = prhs[3];
        const mxArray  *mxHysteresis     = prhs[4];
#if !RFC_MINIMAL
        const mxArray  *mxResidualMethod = prhs[5];
        const mxArray  *mxEnforceMargin  = prhs[6];
        const mxArray  *mxUseHCM         = prhs[7];
        const mxArray  *mxSpreadDamage   = prhs[8];
#endif /*!RFC_MINIMAL*/        

        rfc_value_t    *buffer           = NULL;
        double         *data             = mxGetPr( mxData );
        size_t          data_len         = mxGetNumberOfElements( mxData );
        unsigned        class_count      = (unsigned)( mxGetScalar( mxClassCount ) + 0.5 );
        double          class_width      = mxGetScalar( mxClassWidth );
        double          class_offset     = mxGetScalar( mxClassOffset );
        double          hysteresis       = mxGetScalar( mxHysteresis );
#if !RFC_MINIMAL
        int             residual_method  = (int)( mxGetScalar( mxResidualMethod ) + 0.5 );
        int             enforce_margin   = (int)mxGetScalar( mxEnforceMargin );
        int             use_hcm          = (int)mxGetScalar( mxUseHCM );
        int             spread_damage    = (int)mxGetScalar( mxSpreadDamage );
#else /*RFC_MINIMAL*/
        int             residual_method  = RFC_RES_NONE;
#endif /*!RFC_MINIMAL*/
        size_t          i;
        bool            ok;

        mxAssert( residual_method >= RFC_RES_NONE && residual_method <= RFC_RES_COUNT, 
                  "Invalid residual method!" );

#if !RFC_MINIMAL
#if RFC_HCM_SUPPORT
        mxAssert( use_hcm == 0 || use_hcm == 1,
                  "Invalid HCM flag, use 0 or 1!" );
#else /*!RFC_HCM_SUPPORT*/
        mxAssert( use_hcm == 0,
                  "HCM not supported!" );
#endif /*RFC_HCM_SUPPORT*/


#if RFC_DH_SUPPORT
        mxAssert( spread_damage >= RFC_SD_NONE && spread_damage <= RFC_SD_COUNT,
                  "Invalid spread damage method!" );
#else /*!RFC_DH_SUPPORT*/
        if( spread_damage != 0 )
        {
            mexErrMsgTxt( "Invalid spread damage method, only 0 accepted!" );
        }
#endif /*RFC_DH_SUPPORT*/
#endif /*!RFC_MINIMAL*/

        ok = RFC_init( &rfc_ctx, 
                       class_count, (rfc_value_t)class_width, (rfc_value_t)class_offset, 
                       (rfc_value_t)hysteresis, RFC_FLAGS_DEFAULT );
#if RFC_DEBUG_FLAGS
        rfc_ctx.debug_vfprintf_fcn = rfc_vfprintf_fcn;
        RFC_set_flags( &rfc_ctx, RFC_FLAGS_LOG_CLOSED_CYCLES, 1 );
#endif /*RFC_DEBUG_FLAGS*/

#if RFC_TP_SUPPORT
        if( ok )
        {
            rfc_ctx.tp_cap = 128;
            rfc_ctx.tp     = (rfc_value_tuple_s*)mem_alloc( NULL, rfc_ctx.tp_cap, 
                                                            sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TP );

            ok = RFC_tp_init( &rfc_ctx, rfc_ctx.tp, rfc_ctx.tp_cap, /* is_static */ true );
        }
#endif /*RFC_TP_SUPPORT*/

        if( !ok )
        {
            mexErrMsgTxt( "Error during initialization!" );
        }

#if RFC_DH_SUPPORT
        rfc_ctx.spread_damage_method = spread_damage;
#endif /*RFC_DH_SUPPORT*/

        /* Cast values from double type to rfc_value_t */ 
        if( sizeof( rfc_value_t ) != sizeof(double) && data_len )  /* maybe unsafe! */
        {
            buffer = (rfc_value_t *)mem_alloc( NULL, data_len, 
                                               sizeof(rfc_value_t), RFC_MEM_AIM_TEMP );

            if( !buffer )
            {
                RFC_deinit( &rfc_ctx );
                mexErrMsgTxt( "Error during initialization!" );
            }

            for( i = 0; i < data_len; i++ )
            {
                buffer[i] = (rfc_value_t)data[i];
            }
        }
        else buffer = (rfc_value_t*)data;

        /* Rainflow counting */

#if !RFC_MINIMAL
        /* Setup */
        rfc_ctx.internal.flags  |= enforce_margin ? RFC_FLAGS_ENFORCE_MARGIN : 0;
#endif /*!RFC_MINIMAL*/
#if RFC_HCM_SUPPORT
        rfc_ctx.counting_method  = use_hcm ? RFC_COUNTING_METHOD_HCM : RFC_COUNTING_METHOD_4PTM;
#else /*!RFC_HCM_SUPPORT*/
#if !RFC_MINIMAL
        rfc_ctx.counting_method  = RFC_COUNTING_METHOD_4PTM;
#endif /*!RFC_MINIMAL*/
#endif /*RFC_HCM_SUPPORT*/
        ok = RFC_feed( &rfc_ctx, buffer, data_len ) &&
             RFC_finalize( &rfc_ctx, residual_method );

        /* Free temporary buffer (cast) */
        if( (void*)buffer != (void*)data )
        {
            mem_alloc( buffer, 0, 0, RFC_MEM_AIM_TEMP );
            buffer = NULL;
        }

        if( !ok )
        {
            int error = rfc_ctx.error;

            RFC_deinit( &rfc_ctx );
            switch( error )
            {
                case RFC_ERROR_INVARG:
                    mexErrMsgTxt( "Invalid argument(s)!" );
                case RFC_ERROR_MEMORY:
                    mexErrMsgTxt( "Error during memory allocation!" );
#if RFC_AT_SUPPORT
                case RFC_ERROR_AT:
                    mexErrMsgTxt( "Error during amplitude transformation!" );
#endif /*RFC_AT_SUPPORT*/
#if RFC_TP_SUPPORT
                case RFC_ERROR_TP:
                    mexErrMsgTxt( "Error during turning point access!" );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DAMAGE_FAST
                case RFC_ERROR_LUT:
                    mexErrMsgTxt( "Error during lookup table access!" );
#endif /*RFC_DAMAGE_FAST*/
                case RFC_ERROR_UNEXP:
                default:
                    mexErrMsgTxt( "Unexpected error occured!" );
            }
        }

        /* Return results */
        if( plhs )
        {
            /* Damage */
            plhs[0] = mxCreateDoubleScalar( rfc_ctx.damage );

            /* Residue */
#if RFC_HCM_SUPPORT
            if( use_hcm )
            {
                if( nlhs > 1 && rfc_ctx.internal.hcm.stack )
                {
                    mxArray* re = mxCreateDoubleMatrix( rfc_ctx.internal.hcm.IZ, 1, mxREAL );
                    if( re )
                    {
                        size_t i;
                        double *val = mxGetPr(re);

                        for( i = 0; i < rfc_ctx.internal.hcm.IZ; i++ )
                        {
                            *val++ = (double)rfc_ctx.internal.hcm.stack[i].value;
                        }
                        plhs[1] = re;
                    }
                }
            }
            else
#endif /*!RFC_HCM_SUPPORT*/
            {
                if( nlhs > 1 && rfc_ctx.residue )
                {
                    mxArray* re = mxCreateDoubleMatrix( rfc_ctx.residue_cnt, 1, mxREAL );
                    if( re )
                    {
                        size_t i;
                        double *val = mxGetPr(re);

                        for( i = 0; i < rfc_ctx.residue_cnt; i++ )
                        {
                            *val++ = (double)rfc_ctx.residue[i].value;
                        }
                        plhs[1] = re;
                    }
                }
            }

            /* Rainflow matrix (column major order) */
            if( nlhs > 2 && rfc_ctx.rfm )
            {
                mxArray* rfm = mxCreateDoubleMatrix( class_count, class_count, mxREAL );
                if( rfm )
                {
                    double *ptr = mxGetPr(rfm);
                    size_t from, to;
                    for( to = 0; to < class_count; to++ )
                    {
                        for( from = 0; from < class_count; from++ )
                        {
                            *ptr++ = (double)rfc_ctx.rfm[ from * class_count + to ] / rfc_ctx.full_inc;
                        }
                    }
                    plhs[2] = rfm;
                }
            }
            
#if !RFC_MINIMAL
            /* Range pair */
            if( nlhs > 3 && rfc_ctx.rp )
            {
                mxArray* rp = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( rp )
                {
                    double *ptr = mxGetPr(rp);
                    size_t i;
                    for( i = 0; i < class_count; i++ )
                    {
                        *ptr++ = (double)rfc_ctx.rp[i];
                    }
                    plhs[3] = rp;
                }
            }

            /* Level crossing */
            if( nlhs > 4 && rfc_ctx.lc )
            {
                mxArray* lc = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( lc )
                {
                    double *ptr = mxGetPr(lc);
                    size_t i;
                    for( i = 0; i < class_count; i++ )
                    {
                        *ptr++ = (double)rfc_ctx.lc[i];
                    }
                    plhs[4] = lc;
                }
            }
#if RFC_TP_SUPPORT
            /* Turning points */
            if( nlhs > 5 && rfc_ctx.tp )
            {
#if RFC_DH_SUPPORT
                mxArray *tp  = mxCreateDoubleMatrix( rfc_ctx.tp_cnt, 3, mxREAL );
                double  *dam = tp ? ( mxGetPr(tp) + 2 * rfc_ctx.tp_cnt ) : NULL;
#else /*!RFC_DH_SUPPORT*/
                mxArray* tp  = mxCreateDoubleMatrix( rfc_ctx.tp_cnt, 2, mxREAL );
#endif /*RFC_DH_SUPPORT*/

                if( tp )
                {
                    size_t  i;
                    double *idx  = mxGetPr(tp) + 0 * rfc_ctx.tp_cnt;
                    double *val  = mxGetPr(tp) + 1 * rfc_ctx.tp_cnt;
                    double  D    = 0.0;

                    for( i = 0; i < rfc_ctx.tp_cnt; i++ )
                    {
                        *val++  = (double)rfc_ctx.tp[i].value;
                        *idx++  = (double)rfc_ctx.tp[i].pos;
#if RFC_DH_SUPPORT
                        *dam++  = (double)rfc_ctx.tp[i].damage;
                         D     += (double)rfc_ctx.tp[i].damage;
#endif /*RFC_DH_SUPPORT*/
                    }
                    //assert( D == rfc_ctx.damage );
                    plhs[5] = tp;
                }
            }
#endif /*RFC_TP_SUPPORT*/
#endif /*!RFC_MINIMAL*/
        }

        /* Deinitialize rainflow context */
        RFC_deinit( &rfc_ctx );
    }
}


#if RFC_AT_SUPPORT
/**
 * MATLAB wrapper for the amplitude transformation
 */
static
void mexAmpTransform( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    if( nrhs != 5 )
    {
        mexErrMsgTxt( "Function needs exact 5 arguments!" );
    }
    else
    {
        mwSize n, cnt;
        mxArray *mxResult = NULL;

        rfc_ctx_s ctx = { sizeof(ctx) };

        const mxArray *mxSa          = prhs[0];
        const mxArray *mxSm          = prhs[1];
        const mxArray *mxM           = prhs[2];
        const mxArray *mxTarget      = prhs[3];
        const mxArray *mxTarget_is_R = prhs[4];

        double M           =      mxGetScalar( mxM );
        double target      =      mxGetScalar( mxTarget );
        bool   target_is_R = (int)mxGetScalar( mxTarget_is_R );

        cnt = mxGetNumberOfElements( mxSa );
        mxAssert( mxGetNumberOfElements( mxSm ) == cnt, "Sa and Sm must have same length!" );

        mxResult = mxCreateDoubleMatrix( mxGetDimensions( mxSa )[0], mxGetDimensions( mxSa )[1], mxREAL );
        mxAssert( mxResult, "Memory error!" );

        mxAssert( RFC_init( &ctx, 0 /*class_count*/, 0.0 /*class_width*/, 0.0 /*class_offset*/, 0.0 /*hysteresis*/, RFC_FLAGS_DEFAULT ),
                  "RFC initialization error!" );
        mxAssert( RFC_at_init( &ctx, NULL /*Sa*/, NULL /*Sm*/, 0 /*count*/, M, target /*Sm_rig*/, target /*R_rig*/, target_is_R, false /*symmetric*/ ), 
                  "RFC initialization error!" );

        for( n = 0; n < cnt; n++ )
        {
            double Sa_n = mxGetPr( mxSa )[n];
            double Sm_n = mxGetPr( mxSm )[n];
            double Sa_t;

            RFC_at_transform( &ctx, Sa_n, Sm_n, &Sa_t );
            mxGetPr( mxResult )[n] = Sa_t;
        }

        plhs[0] = mxResult;
    }
}
#endif /*RFC_AT_SUPPORT*/


#if RFC_TP_SUPPORT
/**
 * MATLAB wrapper calculates turning points from data points
 */
static
void mexTP( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    if( nrhs != 3 )
    {
        mexErrMsgTxt( "Function needs exact 3 arguments!" );
    }
    else
    {
        rfc_value_tuple_s *tp = NULL;
        mxArray *mxResult = NULL;
        double *ptr;
        rfc_ctx_s ctx = { sizeof(ctx) };
        mwSize n;
        bool ok = true;

        const mxArray *mxData          = prhs[0];
        const mxArray *mxHysteresis    = prhs[1];
        const mxArray *mxEnforceMargin = prhs[2];

        double hysteresis     = mxGetScalar( mxHysteresis );
        bool   enforce_margin = (int)mxGetScalar( mxEnforceMargin );

        tp = (rfc_value_tuple_s*)calloc( mxGetNumberOfElements( mxData ), sizeof(rfc_value_tuple_s) );
        if( !tp )
        {
            mexErrMsgTxt( "Memory allocation error!" );
        }

        if( !RFC_init( &ctx, 0 /*class_count*/, 0.0 /*class_width*/, 0.0 /*class_offset*/, hysteresis, RFC_FLAGS_DEFAULT ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC init!" );
        }

        if( !RFC_tp_init( &ctx, tp, mxGetNumberOfElements( mxData ), true /*tp_is_static*/ ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC tp init!" );
        }

        if( !RFC_feed( &ctx, mxGetPr( mxData ), mxGetNumberOfElements( mxData ) ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC feed!" );
        }

        if( !finalize_res_ignore( &ctx, ctx.internal.flags ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC finalize!" );
        }

        mxResult = mxCreateDoubleMatrix( 2, ctx.tp_cnt, mxREAL );
        mxAssert( mxResult, "Memory allocation error!" );
        for( ptr = mxGetPr( mxResult ), n = 0; n < ctx.tp_cnt; n++ )
        {
            *ptr++ = tp[n].value;
            *ptr++ = (double)tp[n].pos;
        }
        free( tp );

        if( !RFC_deinit( &ctx ) )
        {
            mexErrMsgTxt( "Error on RFC deinit!" );
        }

        plhs[0] = mxResult;
    }
}
#endif /*RFC_TP_SUPPORT*/


#if !RFC_MINIMAL
/**
 * @brief      Compare two string case insensitive
 *
 * @param[in]  a     First string
 * @param[in]  b     Second string
 *
 * @return     0 on equality
 */
static
int wal_stricmp( const char *a, const char *b )
{
    int ca, cb;
    do
    {
        ca = (unsigned char) *a++;
        cb = (unsigned char) *b++;
        ca = tolower( toupper(ca) );
        cb = tolower( toupper(cb) );
    }
    while( ca == cb && ca != '\0' );
    
    return ca - cb;
}


/**
 * The MATLAB MEX main function
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    char buffer[80];

    if( !nrhs || !mxIsChar( prhs[0] ) || 0 != mxGetString( prhs[0], buffer, sizeof(buffer) ) )
    {
        mexPrintf( "%s", RFC_MEX_USAGE );
        return;
    }

    if( 0 == wal_stricmp( buffer, "rfc" ) )
    { 
        mexRainflow( nlhs, plhs, nrhs - 1, prhs + 1 );
    }
#if RFC_AT_SUPPORT
    else if( 0 == wal_stricmp( buffer, "amptransform" ) )
    {
        mexAmpTransform( nlhs, plhs, nrhs - 1, prhs + 1 );
    }
#endif /*RFC_AT_SUPPORT*/
#if RFC_TP_SUPPORT
    else if( 0 == wal_stricmp( buffer, "turningpoints" ) )
    {
        mexTP( nlhs, plhs, nrhs - 1, prhs + 1 );
    }
#endif /*RFC_TP_SUPPORT*/
    else
    {
        mexPrintf( "Unknown subfunction \"%s\"!\n", buffer );
    }
}
#else /*RFC_MINIMAL*/
/**
 * The MATLAB MEX main function
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    mexRainflow( nlhs, plhs, nrhs, prhs );
}
#endif /*!RFC_MINIMAL*/

#endif /*MATLAB_MEX_FILE*/
