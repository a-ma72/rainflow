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


#include "rainflow.h"

#include <assert.h>  /* assert() */
#include <math.h>    /* exp(), log(), fabs() */
#include <string.h>  /* memset() */


#ifdef MATLAB_MEX_FILE
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
static void                 RFC_reset                           ( rfc_ctx_s * );
static void                 RFC_cycle_find                      ( rfc_ctx_s *, int flags );
#else /*RFC_MINIMAL*/
#define mexFunction         mexRainflow
#define RFC_cycle_find      RFC_cycle_find_4ptm
#endif /*!RFC_MINIMAL*/
static bool                 RFC_feed_once                       ( rfc_ctx_s *, const rfc_value_tuple_s* tp, int flags );
#if RFC_DH_SUPPORT
static bool                 RFC_feed_once_dh                    ( rfc_ctx_s * );
#endif /*RFC_DH_SUPPORT*/
#if RFC_TP_SUPPORT
static bool                 RFC_feed_once_tp_check_margin       ( rfc_ctx_s *, const rfc_value_tuple_s* pt, rfc_value_tuple_s** tp_residue );
#endif /*RFC_TP_SUPPORT*/
static bool                 RFC_feed_finalize                   ( rfc_ctx_s * );
#if RFC_TP_SUPPORT
static bool                 RFC_feed_finalize_tp                ( rfc_ctx_s *, rfc_value_tuple_s *tp_interim, int flags );
#endif /*RFC_TP_SUPPORT*/
#if RFC_HCM_SUPPORT
static bool                 RFC_feed_finalize_hcm               ( rfc_ctx_s *, int flags );
#endif /*!RFC_HCM_SUPPORT*/
static rfc_value_tuple_s *  RFC_feed_filter_pt                  ( rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 RFC_cycle_find_4ptm                 ( rfc_ctx_s *, int flags );
#if RFC_HCM_SUPPORT
static void                 RFC_cycle_find_hcm                  ( rfc_ctx_s *, int flags );
#endif /*RFC_HCM_SUPPORT*/
#if !RFC_MINIMAL
static void                 RFC_cycle_process_lc                ( rfc_ctx_s *, int flags );
#endif /*!RFC_MINIMAL*/
static void                 RFC_cycle_process_counts            ( rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, int flags );
/* Methods on residue */
static bool                 RFC_finalize_res_ignore             ( rfc_ctx_s *, int flags );
#if !RFC_MINIMAL
static bool                 RFC_finalize_res_discard            ( rfc_ctx_s *, int flags );
static bool                 RFC_finalize_res_weight_cycles      ( rfc_ctx_s *, RFC_counts_type weight, int flags );
static bool                 RFC_finalize_res_clormann_seeger    ( rfc_ctx_s *, int flags );
static bool                 RFC_finalize_res_rp_DIN45667        ( rfc_ctx_s *, int flags );
static bool                 RFC_finalize_res_repeated           ( rfc_ctx_s *, int flags );
static bool                 RFC_residue_exchange                ( rfc_ctx_s *, rfc_value_tuple_s **residue, size_t *residue_cap, size_t *residue_cnt, bool restore );
#endif /*!RFC_MINIMAL*/
static void                 RFC_residue_remove_item             ( rfc_ctx_s *, size_t index, size_t count );
/* Memory allocator */
static void *               RFC_mem_alloc                       ( void *ptr, size_t num, size_t size, int aim );
#if RFC_TP_SUPPORT
/* Methods on turning points history */
static bool                 RFC_tp_add                          ( rfc_ctx_s *, rfc_value_tuple_s *pt );
static void                 RFC_tp_lock                         ( rfc_ctx_s *, bool do_lock );
static void                 RFC_tp_refeed                       ( rfc_ctx_s *, RFC_value_type new_hysteresis, const rfc_class_param_s *new_class_param );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
static void                 RFC_spread_damage                   ( rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, int flags );
#endif /*RFC_DH_SUPPORT*/
#if RFC_AT_SUPPORT
static double               RFC_at_R_to_Sm_norm                 ( rfc_ctx_s *, double R );
static double               RFC_at_alleviation                  ( rfc_ctx_s *, double Sm_norm );
#endif /*RFC_AT_SUPPORT*/
/* Other */
#if !RFC_MINIMAL
static void                 RFC_class_param_set                 ( rfc_ctx_s *, const rfc_class_param_s * );
static void                 RFC_class_param_get                 ( rfc_ctx_s *, rfc_class_param_s * );
#endif /*!RFC_MINIMAL*/
static bool                 RFC_error_raise                     ( rfc_ctx_s *, int );
static double               RFC_damage_calc_amplitude           ( rfc_ctx_s *, double amplitude );
static double               RFC_damage_calc                     ( rfc_ctx_s *, unsigned class_from, unsigned class_to );
#if RFC_DAMAGE_FAST
static void                 RFC_damage_lut_init                 ( rfc_ctx_s * );
static double               RFC_damage_calc_fast                ( rfc_ctx_s *, unsigned class_from, unsigned class_to );
#endif /*RFC_DAMAGE_FAST*/
static RFC_value_type       value_delta                         ( RFC_value_type from, RFC_value_type to, int *sign_ptr );


#define QUANTIZE( r, v )   ( (unsigned)( ((v) - (r)->class_offset) / (r)->class_width ) )
#define CLASS_MEAN( r, c ) ( (double)( (r)->class_width * (0.5 + (c)) + (r)->class_offset ) )
#define NUMEL( x )         ( sizeof(x) / sizeof(*(x)) )


/**
 * @brief      Initialization (rainflow context).
 *
 * @param      ctx           The rainflow context
 * @param[in]  class_count   The class count
 * @param[in]  class_width   The class width
 * @param[in]  class_offset  The class offset
 * @param[in]  hysteresis    The hysteresis
#if RFC_TP_SUPPORT 
 * @param[in]  tp            Pointer to turning points buffer
 * @param[in]  tp_cap        Number of turning points in buffer
#endif
 *
 * @return     true on success
 */
bool RFC_init                 ( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                                           RFC_value_type hysteresis )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;
        return false;
    }

    if( rfc_ctx->state != RFC_STATE_INIT0 )
    {
        return false;
    }

    /* Flags */
#if !RFC_MINIMAL
    rfc_ctx->flags                      = RFC_FLAGS_COUNT_ALL            | 
#if RFC_TP_SUPPORT
                                          RFC_FLAGS_TPPRUNE_PRESERVE_POS | 
                                          RFC_FLAGS_TPPRUNE_PRESERVE_RES |
#endif /*RFC_TP_SUPPORT*/
                                          0;
#else /*RFC_MINIMAL*/
    rfc_ctx->flags                      = RFC_FLAGS_COUNT_MATRIX | 
                                          RFC_FLAGS_COUNT_DAMAGE;
#endif /*!RFC_MINIMAL*/
    RFC_flags_set( rfc_ctx, rfc_ctx->flags );  /* Validate flags */

    /* Counter increments */
    rfc_ctx->full_inc                   = RFC_FULL_CYCLE_INCREMENT;
    rfc_ctx->half_inc                   = RFC_HALF_CYCLE_INCREMENT;
    rfc_ctx->curr_inc                   = RFC_FULL_CYCLE_INCREMENT;

    if( class_count )
    {
        if( class_count > 512 || class_width <= 0.0 )
        {
            rfc_ctx->error = RFC_ERROR_INVARG;
            return false;
        }
    }

    /* Rainflow class parameters */
    rfc_ctx->class_count                = class_count;
    rfc_ctx->class_width                = class_width;
    rfc_ctx->class_offset               = class_offset;
    rfc_ctx->hysteresis                 = hysteresis;

    /* Woehler curve (fictive) */
    rfc_ctx->wl_sd                      =  1e3;            /* Fictive amplitude */
    rfc_ctx->wl_nd                      =  1e7;            /* Fictive count */
    rfc_ctx->wl_k                       = -5.0;            /* Fictive gradient */
#if !RFC_MINIMAL
    rfc_ctx->wl_k2                      =  rfc_ctx->wl_k;  /* "Miner elementar", if k == k2 */
    rfc_ctx->wl_omission                =  0.0;            /* No omission per default */
#endif /*!RFC_MINIMAL*/

    /* Memory allocator */
    if( !rfc_ctx->mem_alloc )
    {
        rfc_ctx->mem_alloc = RFC_mem_alloc;
    }
    
#if RFC_USE_DELEGATES
    /* Delegates (optional, set to NULL for standard or to your own functions! ) */
#if RFC_TP_SUPPORT
    rfc_ctx->tp_next_fcn                    = NULL;
    rfc_ctx->tp_add_fcn                     = NULL;
#endif /*RFC_TP_SUPPORT*/
    rfc_ctx->cycle_find_fcn                 = NULL;
    rfc_ctx->finalize_fcn                   = NULL;
    rfc_ctx->damage_calc_fcn                = NULL;
#endif /*RFC_USE_DELEGATES*/

#if !RFC_MINIMAL
    /* Rainflow counting method */
    rfc_ctx->counting_method                = RFC_COUNTING_METHOD_4PTM;
#endif /*RFC_MINIMAL*/

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
        /* Non-sparse storages (optional, may be NULL) */
        rfc_ctx->matrix                     = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count * class_count, 
                                                                                    sizeof(RFC_counts_type), RFC_MEM_AIM_MATRIX );
#if !RFC_MINIMAL
        rfc_ctx->rp                         = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count,
                                                                                    sizeof(RFC_counts_type), RFC_MEM_AIM_RP );
        rfc_ctx->lc                         = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count,
                                                                                    sizeof(RFC_counts_type), RFC_MEM_AIM_LC );
#endif /*!RFC_MINIMAL*/
    }

    /* Damage */
    rfc_ctx->pseudo_damage              = 0.0;

    /* Internals */
    rfc_ctx->internal.slope             = 0;
    rfc_ctx->internal.extrema[0]        = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]        = nil;  /* local maximum */
#if RFC_GLOBAL_EXTREMA    
    rfc_ctx->internal.extrema_changed   = false;
#endif /*RFC_GLOBAL_EXTREMA*/
#if RFC_TP_SUPPORT
    rfc_ctx->internal.margin[0]         = nil;  /* left  margin */
    rfc_ctx->internal.margin[1]         = nil;  /* right margin */
    rfc_ctx->internal.margin_stage      = 0;
#endif /*RFC_TP_SUPPORT*/

#if !RFC_MINIMAL
    if( rfc_ctx->class_count && ( !rfc_ctx->residue || !rfc_ctx->matrix || !rfc_ctx->rp || !rfc_ctx->lc ) )
#else /*RFC_MINIMAL*/
    if( rfc_ctx->class_count && ( !rfc_ctx->residue || !rfc_ctx->matrix ) )
#endif /*!RFC_MINIMAL*/
    {
        RFC_deinit( rfc_ctx );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }
    
#if RFC_TP_SUPPORT
    /* Turning points storage (optional, may be NULL) */
    rfc_ctx->tp                         = NULL;
    rfc_ctx->tp_cap                     = 0;
    rfc_ctx->tp_cnt                     = 0;
    rfc_ctx->tp_locked                  = 0;
    rfc_ctx->tp_prune_threshold         = (size_t)-1;
    rfc_ctx->tp_prune_size              = (size_t)-1;
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
    rfc_ctx->at.Sa                      = NULL;
    rfc_ctx->at.Sm                      = NULL;
    rfc_ctx->at.count                   = 0;
    rfc_ctx->at.M                       = 0.0;
    rfc_ctx->at.Sm_rig                  = 0.0;
    rfc_ctx->at.R_rig                   = 0.0;
    rfc_ctx->at.R_pinned                = false;

    rfc_ctx->internal.at.count          = 0;
#endif /*RFC_AT_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT;

#if RFC_DAMAGE_FAST
    if( rfc_ctx->class_count )
    {
        rfc_ctx->damage_lut                 = (double*)rfc_ctx->mem_alloc( rfc_ctx->damage_lut, class_count * class_count, 
                                                                           sizeof(double), RFC_MEM_AIM_DLUT );
        RFC_damage_lut_init( rfc_ctx );
    }
#endif /*RFC_DAMAGE_FAST*/

    return true;
}

#if RFC_TP_SUPPORT
/**
 * @brief      Initialize tp buffer
 *
 * @param      ctx        The rainflow context
 * @param[out] tp         The buffer for tp
 * @param[in]  tp_cap     The tp capability
 * @param[in]  is_static  Indicates if tp is static
 *
 * @return     true on success
 */
bool RFC_tp_init( void *ctx, rfc_value_tuple_s *tp, size_t tp_cap, bool is_static )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) || !tp )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
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


bool RFC_tp_init_autoprune( void *ctx, bool autoprune, size_t size, size_t threshold )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->flags                      = rfc_ctx->flags & ~RFC_FLAGS_TPAUTOPRUNE | (autoprune ? RFC_FLAGS_TPAUTOPRUNE : 0);
    rfc_ctx->tp_prune_threshold         = threshold;
    rfc_ctx->tp_prune_size              = size;
}


/**
 * @brief      Drop turning points from storage, to avoid memory excess
 *
 * @param      ctx    The rainflow context
 * @param[in]  limit  The excepted number of points left in turning points
 *                    storage (May be more, if residuals aren't neglected)
 * @param[in]  flags  The flags (see RFC_FLAGS_TPPRUNE_...)
 *
 * @return     true on success
 */
bool RFC_tp_prune( void *ctx, size_t limit, int flags )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) || !rfc_ctx->tp )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

#if RFC_USE_DELEGATES
    if( rfc_ctx->tp_prune_fcn )
    {
        return rfc_ctx->tp_prune_fcn( rfc_ctx, limit, flags );
    }
#endif /*RFC_USE_DELEGATES*/

    if( rfc_ctx->tp && rfc_ctx->tp_cnt > limit )
    {
        rfc_value_tuple_s   *src_beg_it,    /* Source (begin) */
                            *src_end_it,    /* Source (end) */
                            *src_it,        /* Source iterator */
                            *dst_it,        /* Destination iterator */
                            *res_it;        /* Residue iterator */
        size_t               src_pos,       /* Source, position base 1 */
                             dst_i,         /* New turning points, index base 0 */
                             res_i;         /* Residue, index base 0 */

        size_t               removal;       /* Number of turning points to remove */
        size_t               offset;        /* Position offset (minuend) */
        bool                 preserve_pos;  /* Don't justify position */
        bool                 preserve_res;  /* Don't remove turning points, if referenced by resiude */

        removal     = rfc_ctx->tp_cnt - limit;
        dst_it      = rfc_ctx->tp;
        dst_i       = 0;
        src_beg_it  = dst_it + removal;
        src_end_it  = rfc_ctx->tp + rfc_ctx->tp_cnt
                      + ( ( rfc_ctx->state == RFC_STATE_BUSY_INTERIM ) ? 1 : 0 );
        src_it      = src_beg_it;
        src_pos     = removal + 1;
        res_it      = rfc_ctx->residue;
        res_i       = 0;
        offset      = 0;

        preserve_pos = ( flags & RFC_FLAGS_TPPRUNE_PRESERVE_POS ) > 0;  /* Preserve (stream) position */
        preserve_res = ( flags & RFC_FLAGS_TPPRUNE_PRESERVE_RES ) > 0;  /* Preserve residual turning points */

        /* Move turning points ahead */
        while( src_it < src_end_it || res_i < rfc_ctx->residue_cnt )
        {
            /* Check if there are still residual points to consider */
            while( res_i < rfc_ctx->residue_cnt && res_it->tp_pos <= src_pos )
            {
                /* Check if residue refers a turning point from removal area */

                /* First new turning point delivers new offset */
                if( !res_i && !preserve_pos )
                {
                    offset = res_it->pos;
                    assert( offset );
                    offset--;
                }

                /* Residual point is at current source position? */
                if( res_it->tp_pos == src_pos )
                {
                    src_it++;
                    src_pos++;
                }

                if( preserve_res )
                {
                    /* Adjust residue reference information */
                    res_it->tp_pos = dst_i + 1;
                    res_it->pos   -= offset;
                    *dst_it++      = *res_it++;
                    dst_i++;
                    res_i++;
                }
                else
                {
                    /* Residual turning point refers first point now */
                    res_it->tp_pos = 0;  /* Index 0 => "none" */
                    res_it->pos   -= offset;
                    res_it++;
                    res_i++;
                }
            }

            if( src_it < src_end_it )
            {
                /* First new turning point delivers new offset */
                if( !dst_i && !preserve_pos )
                {
                    offset = src_it->pos;
                    assert( offset );
                    offset--;
                }

                /* Copy turning point from source */
                src_it->tp_pos = dst_i + 1;
                src_it->pos   -= offset;
                *dst_it++      = *src_it++;
                dst_i++;
                src_pos++;
            }
        }

        rfc_ctx->tp_cnt                  = dst_i;
        rfc_ctx->internal.pos           -= offset;
        rfc_ctx->internal.global_offset += offset;
    }
    
    return true;
}
#endif /*RFC_TP_SUPPORT*/


#if RFC_DH_SUPPORT
/**
 * @brief      Initialize damage history storage
 *
 * @param      ctx        The rainflow context
 * @param[out] dh         The storage buffer
 * @param[in]  dh_cap     The capacity of dh
 * @param[in]  is_static  true, if dh is static and should not be freed
 *
 * @return     true on success
 */
bool RFC_dh_init( void *ctx, double *dh, size_t dh_cap, bool is_static )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) || rfc_ctx->dh )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state != RFC_STATE_INIT )
    {
        return false;
    }

    rfc_ctx->dh     = dh;
    rfc_ctx->dh_cap = dh_cap;
    rfc_ctx->dh_cnt = 0;

    rfc_ctx->internal.dh_static = is_static;

    return true;
}
#endif /*RFC_DH_SUPPORT*/


#if RFC_AT_SUPPORT
/**
 * @brief      Initialize amplitude transformation
 *
 * @param      ctx        The rainflow context
 * @param[in]  Sa         The reference curve vector, amplitude part
 * @param[in]  Sm         The reference curve vector, mean load part. If Sa and
 *                        Sm_norm are NULL, the standard (FKM) is applied
 * @param[in]  count      The capacity of Sa and Sm
 * @param[in]  M          The mean stress sensitivity
 * @param[in]  Sm_rig     The mean load applied on the test rig
 * @param[in]  R_rig      The mean load ratio applied on the test rig
 * @param[in]  R_pinned   true, if R is constant on test rig (R_rig is used).
 *                        false if Sm is constant on test rig (Sm_rig is used)
 * @param[in]  symmetric  true if Haigh diagram is symmetric at Sa(R=-1)
 *
 * @return     true on success
 */
bool RFC_at_init( void *ctx, const double *Sa, const double *Sm, unsigned count, 
                             double M, double Sm_rig, double R_rig, bool R_pinned, bool symmetric )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) || M < 0.0 )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
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
            rfc_ctx->error = RFC_ERROR_INVARG;
            return false;
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
            rfc_ctx->error = RFC_ERROR_INVARG;
            return false;
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

        assert( !Sa && !Sm && !count );

        Sa_R_Inf = 1.0 / ( 1.0 - M );                      /* y = -x && y = Sa(R=-1) - Mx                  */
        Sa_R_0   = 1.0 / ( 1.0 + M );                      /* y =  x && y = Sa(R=-1) - Mx                  */
        Sa_R_0p5 = Sa_R_0 * ( 1.0 + M/3 ) / ( 1.0 + M );   /* Backtrace Sa(R=0) to Sa(R=-1) with M/3, then */
                                                           /* 3y = x && y = Sa(R=-1) - (M/3)x              */
        if( symmetric )
        {
            /* Build symmetrical reference curve */
            /* Symmetric around R=-1 (Sm=0) */

            double *Sa_ = rfc_ctx->internal.at.Sa;
            double *Sm_ = rfc_ctx->internal.at.Sm;

            assert( NUMEL( rfc_ctx->internal.at.Sa ) >= 5 );
            
            rfc_ctx->internal.at.count = 5;

            Sa_[0] = Sa_R_0p5; Sm_[0] = -Sa_R_0p5 * 3.0;
            Sa_[1] = Sa_R_0;   Sm_[1] = -Sa_R_0;
            Sa_[2] = 1.0;      Sm_[2] =  0.0;
            Sa_[3] = Sa_[1];   Sm_[3] = -Sm_[1];
            Sa_[4] = Sa_[0];   Sm_[4] = -Sm_[0];
        }
        else
        {
            /* Build non-symmetric reference curve */
            double *Sa_ = rfc_ctx->internal.at.Sa;
            double *Sm_ = rfc_ctx->internal.at.Sm;

            assert( NUMEL( rfc_ctx->internal.at.Sa ) >= 3 );
            
            rfc_ctx->internal.at.count = 3;

            Sa_[0] = Sa_R_Inf; Sm_[0] = -Sa_R_Inf;
            Sa_[1] = Sa_R_0;   Sm_[1] =  Sa_R_0;
            Sa_[2] = Sa_R_0p5; Sm_[2] =  Sa_R_0p5 * 3.0;
        }

        rfc_ctx->at.Sa       = rfc_ctx->internal.at.Sa;
        rfc_ctx->at.Sm       = rfc_ctx->internal.at.Sm;
        rfc_ctx->at.count    = rfc_ctx->internal.at.count;
        rfc_ctx->at.M        = M;
        rfc_ctx->at.Sm_rig   = Sm_rig;
        rfc_ctx->at.R_rig    = R_rig;
        rfc_ctx->at.R_pinned = R_pinned;
    }

#if RFC_DAMAGE_FAST
    RFC_damage_lut_init( rfc_ctx );
#endif /*RFC_DAMAGE_FAST*/

    return true;
}
#endif /*RFC_AT_SUPPORT*/


/**
 * @brief      De-initialization (freeing memory).
 *
 * @param      ctx   The rainflow context
 *
 * @return     true on success
 */
bool RFC_deinit( void *ctx )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    if( !rfc_ctx->internal.res_static &&
        rfc_ctx->residue )              rfc_ctx->mem_alloc( rfc_ctx->residue,    0, 0, RFC_MEM_AIM_RESIDUE );
    if( rfc_ctx->matrix )               rfc_ctx->mem_alloc( rfc_ctx->matrix,     0, 0, RFC_MEM_AIM_MATRIX );
#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )           rfc_ctx->mem_alloc( rfc_ctx->damage_lut, 0, 0, RFC_MEM_AIM_DLUT );
#endif /*RFC_DAMAGE_FAST*/
#if !RFC_MINIMAL
    if( rfc_ctx->rp )                   rfc_ctx->mem_alloc( rfc_ctx->rp,         0, 0, RFC_MEM_AIM_RP );
    if( rfc_ctx->lc )                   rfc_ctx->mem_alloc( rfc_ctx->lc,         0, 0, RFC_MEM_AIM_LC );
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
    if( rfc_ctx->tp && !rfc_ctx->internal.tp_static )
    {
                                        rfc_ctx->mem_alloc( rfc_ctx->tp,         0, 0, RFC_MEM_AIM_TP );
    }           
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
    if( rfc_ctx->dh && !rfc_ctx->internal.dh_static )
    {               
                                        rfc_ctx->mem_alloc( rfc_ctx->dh,         0, 0, RFC_MEM_AIM_DH );
    }
#endif /*RFC_DH_SUPPORT*/

#if RFC_DAMAGE_FAST
    rfc_ctx->damage_lut                 = NULL;
#endif /*RFC_DAMAGE_FAST*/

    rfc_ctx->residue                    = NULL;
    rfc_ctx->residue_cap                = 0;
    rfc_ctx->residue_cnt                = 0;

    rfc_ctx->matrix                     = NULL;
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
    rfc_ctx->internal.global_offset     = 0;
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

    rfc_ctx->internal.at.count          = 0;
#endif /*RFC_AT_SUPPORT*/

#if RFC_HCM_SUPPORT
    /* Remove stack */
    if( rfc_ctx->internal.hcm.stack )   rfc_ctx->mem_alloc( rfc_ctx->internal.hcm.stack, 0, 0, RFC_MEM_AIM_HCM );

    rfc_ctx->internal.hcm.stack         = NULL;
    rfc_ctx->internal.hcm.stack_cap     = 0;

    /* Stack pointers */
    rfc_ctx->internal.hcm.IZ         = 0;
    rfc_ctx->internal.hcm.IR         = 1;

#endif /*RFC_HCM_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT0;

    return true;
}


/**
 * @brief      Set counting flags
 *
 * @param      ctx    The rainflow context
 * @param[in]  flags  The flags
 *
 * @return     Old flags
 */
static
int RFC_flags_set( void *ctx, int flags )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;
    int old_flags;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    if( !rfc_ctx->class_count )
    {
        flags &= ~(RFC_FLAGS_COUNT_ALL);
    }

    old_flags = rfc_ctx->flags;
    rfc_ctx->flags = flags;

    return old_flags;
}


/**
 * @brief      "Feed" counting algorithm with data samples (consecutive calls
 *             allowed).
 *
 * @param      ctx         The rainflow context
 * @param[in]  data        The data
 * @param[in]  data_count  The data count
 *
 * @return     true on success
 */
bool RFC_feed( void *ctx, const RFC_value_type * data, size_t data_count )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( data_count && !data ) return false;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    /* Process data */
    while( data_count-- )
    {
        rfc_value_tuple_s tp = { *data++ };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

        /* Assign class and global position (base 1) */
        tp.class = QUANTIZE( rfc_ctx, tp.value );
        tp.pos   = ++rfc_ctx->internal.pos;

        if( !RFC_feed_once( rfc_ctx, &tp, rfc_ctx->flags ) ) return false;
    }

    return true;
}


#if !RFC_MINIMAL
/**
 * @brief      Feed counting algorithm with data tuples.
 *
 * @param      ctx         The rainflow context
 * @param[in]  data        The data tuples
 * @param[in]  data_count  The data count
 *
 * @return     true on success
 */
bool RFC_feed_tuple( void *ctx, rfc_value_tuple_s *data, size_t data_count )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( data_count && !data ) return false;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

    /* Process data */
    while( data_count-- )
    {
        if( !RFC_feed_once( rfc_ctx, data++, rfc_ctx->flags ) ) return false;
    }

    return true;
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Finalize pending counts and turning point storage.
 *
 * @param      ctx              The rainflow context
 * @param[in]  residual_method  The residual method (RFC_RES_...)
 *
 * @return     true on success
 */
bool RFC_finalize( void *ctx, int residual_method )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;
    bool ok;

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }
    
    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        return false;
    }

#if RFC_USE_DELEGATES
    if( rfc_ctx->finalize_fcn )
    {
        ok = rfc_ctx->finalize_fcn( rfc_ctx, residual_method );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        int flags = rfc_ctx->flags;

#if !RFC_MINIMAL
        flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

        switch( residual_method )
        {
            case RFC_RES_NONE:
                /* FALLTHROUGH */
            case RFC_RES_IGNORE:
                ok = RFC_finalize_res_ignore( rfc_ctx, flags );
                break;
#if !RFC_MINIMAL
            case RFC_RES_DISCARD:
                ok = RFC_finalize_res_discard( rfc_ctx, flags );
                break;
            case RFC_RES_HALFCYCLES:
                ok = RFC_finalize_res_weight_cycles( rfc_ctx, rfc_ctx->half_inc, flags );
                break;
            case RFC_RES_FULLCYCLES:
                ok = RFC_finalize_res_weight_cycles( rfc_ctx, rfc_ctx->full_inc, flags );
                break;
            case RFC_RES_CLORMANN_SEEGER:
                ok = RFC_finalize_res_clormann_seeger( rfc_ctx, flags );
                break;
            case RFC_RES_REPEATED:
                ok = RFC_finalize_res_repeated( rfc_ctx, flags );
                break;
            case RFC_RES_RP_DIN45667:
                ok = RFC_finalize_res_rp_DIN45667( rfc_ctx, flags );
                break;
#endif /*!RFC_MINIMAL*/
            default:
                assert( false );
                rfc_ctx->error = RFC_ERROR_INVARG;
                ok = false;
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

    rfc_ctx->state = ok ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
    return ok;
}


#if RFC_AT_SUPPORT
/**
 * @brief      Amplitude transformation to take mean load influence into
 *             account.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  Sa       Amplitude
 * @param[in]  Sm       Mean load
 *
 * @return     Transformed amplitude Sa
 */
double RFC_at_transform( rfc_ctx_s *rfc_ctx, double Sa, double Sm )
{
    double Sa_transform = Sa;

    assert( rfc_ctx );

    /* Amplitude is always positive */
    Sa = fabs( Sa );

#if RFC_USE_DELEGATES
    if( rfc_ctx->at_transform_fcn )
    {
        return rfc_ctx->at_transform_fcn( rfc_ctx, Sa, Sm );
    }
#endif

    if (!rfc_ctx->at.count)
    {
        /* No reference curve given, return original amplitude */
        return Sa;
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

        /* Normalize Sm (Sa=1) */
        Sm_norm_base = Sm / Sa;
        alleviation_base = RFC_at_alleviation( rfc_ctx, Sm_norm_base );

        if( rfc_ctx->at.R_pinned )
        {
            /* Calculate intersection of R slope and M slope */
            Sm_norm_target = RFC_at_R_to_Sm_norm( rfc_ctx, rfc_ctx->at.R_rig );
            Sa_transform   = Sa / alleviation_base * RFC_at_alleviation( rfc_ctx, Sm_norm_target );
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
                        /* Next segment */
                        Sa_rhs = Sa / alleviation_base * RFC_at_alleviation( rfc_ctx, Sm_[n] / Sa_[n] );
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
                    /* First segment */
                    Sa_rhs = Sa / alleviation_base * RFC_at_alleviation( rfc_ctx, Sm_[0] / Sa_[0] );
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

    return Sa_transform;
}
#endif /*RFC_AT_SUPPORT*/










/*** Implementation static functions ***/

#if !RFC_MINIMAL
/**
 * brief       Reset data processing information (empty containers)
 *
 * @param      rfc_ctx  The rainflow context
 */
static
void RFC_reset( rfc_ctx_s *rfc_ctx )
{
    rfc_value_tuple_s  nil = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    assert( rfc_ctx && rfc_ctx->state >= RFC_STATE_INIT );

    if( rfc_ctx->matrix )
    {
        memset( rfc_ctx->matrix, 0, sizeof(RFC_counts_type) * rfc_ctx->class_count * rfc_ctx->class_count );
    }

    if( rfc_ctx->rp )
    {
        memset( rfc_ctx->rp, 0, sizeof(RFC_counts_type) * rfc_ctx->class_count );
    }

    if( rfc_ctx->lc )
    {
        memset( rfc_ctx->lc, 0, sizeof(RFC_counts_type) * rfc_ctx->class_count );
    }

    rfc_ctx->residue_cnt                = 0;

    rfc_ctx->internal.slope             = 0;
    rfc_ctx->internal.extrema[0]        = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]        = nil;  /* local maximum */
#if RFC_GLOBAL_EXTREMA
    rfc_ctx->internal.extrema_changed   = false;
#endif
    rfc_ctx->internal.pos               = 0;
    rfc_ctx->internal.global_offset     = 0;
#if RFC_TP_SUPPORT
    rfc_ctx->internal.margin[0]         = nil;  /* left margin */
    rfc_ctx->internal.margin[1]         = nil;  /* right margin */
    rfc_ctx->internal.margin_stage      = 0;
    rfc_ctx->tp_cnt                     = 0;
    rfc_ctx->tp_locked                  = 0;
#endif /*RFC_TP_SUPPORT*/

#if RFC_DH_SUPPORT
    rfc_ctx->dh_cnt                     = 0;
#endif /*RFC_DH_SUPPORT*/

    rfc_ctx->pseudo_damage              = 0.0;

#if RFC_AT_SUPPORT
    rfc_ctx->at.Sa                      = NULL;
    rfc_ctx->at.Sm                      = NULL;
    rfc_ctx->at.count                   = 0;
    rfc_ctx->at.M                       = 0.0;
    rfc_ctx->at.Sm_rig                  = 0.0;
    rfc_ctx->at.R_rig                   = 0.0;
    rfc_ctx->at.R_pinned                = false;

    rfc_ctx->internal.at.count          = 0;
#endif /*RFC_AT_SUPPORT*/

#if RFC_DAMAGE_FAST
    RFC_damage_lut_init( rfc_ctx );
#endif /*RFC_DAMAGE_FAST*/

#if RFC_HCM_SUPPORT
    /* Reset stack pointers */
    rfc_ctx->internal.hcm.IR            = 1;
    rfc_ctx->internal.hcm.IZ            = 0;
#endif /*RFC_HCM_SUPPORT*/

    rfc_ctx->state = RFC_STATE_INIT;
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Processing one data point. Find turning points and check for
 *             closed cycles.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  pt       The data tuple
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_feed_once( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt, int flags )
{
    rfc_value_tuple_s *tp_residue;  /* Pointer to residue element */

    assert( rfc_ctx && pt );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

#if RFC_DH_SUPPORT
    /* Resize damage history if necessary */
    if( !RFC_feed_once_dh( rfc_ctx ) )
    {
        return false;
    }
#endif /*RFC_DH_SUPPORT*/

    /* Check for next turning point and update residue. tp_residue is NULL, if there is no turning point */
    /* Otherwise tp_residue refers the forelast element in member rfc_ctx->residue */
    tp_residue = RFC_feed_filter_pt( rfc_ctx, pt );

#if RFC_TP_SUPPORT
    /* Check if pt influences margins (tp_residue may be set to NULL then!) */
    if( !RFC_feed_once_tp_check_margin( rfc_ctx, pt, &tp_residue ) )
    {
        return false;
    }
#endif /*RFC_TP_SUPPORT*/

    /* Countings */

    /* Add turning point and check for closed cycles */
    if( tp_residue )
    {
#if RFC_TP_SUPPORT
        /* Add new turning point, position in field tp will be stored in tp_residue */
        if( !RFC_tp_add( rfc_ctx, tp_residue ) )
        {
            return false;
        }
#endif /*RFC_TP_SUPPORT*/

#if !RFC_MINIMAL
        /* New turning point, do LC count */
        RFC_cycle_process_lc( rfc_ctx, flags & (RFC_FLAGS_COUNT_LC | RFC_FLAGS_ENFORCE_MARGIN) );
#endif /*!RFC_MINIMAL*/

        if( rfc_ctx->class_count )
        {
#if !RFC_MINIMAL
            flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

            /* Check for closed cycles and count. Modifies residue! */
            RFC_cycle_find( rfc_ctx, flags );
        }
        else
        {
            if( rfc_ctx->residue_cnt > 1 )
            {
                RFC_residue_remove_item( rfc_ctx, 0, 1 );
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
 *
 * @return     true on success
 */
static
bool RFC_feed_once_dh( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->dh && rfc_ctx->dh_cnt > rfc_ctx->dh_cap )
    {
        size_t new_cap = (size_t)1024 * ( rfc_ctx->dh_cap / 640 + 1 ); /* + 60% + 1024 */

        rfc_ctx->dh = (double*)rfc_ctx->mem_alloc( rfc_ctx->dh, new_cap, 
                                                   sizeof(RFC_value_type), RFC_MEM_AIM_DH );

        if( !rfc_ctx->dh )
        {
            return RFC_error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }

        memset( rfc_ctx->dh + rfc_ctx->dh_cnt, 0, sizeof(RFC_value_type) * ( new_cap - rfc_ctx->dh_cap ) );
        rfc_ctx->dh_cap = new_cap;
        rfc_ctx->dh_cnt++;
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
bool RFC_feed_once_tp_check_margin( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt, rfc_value_tuple_s** tp_residue )
{
    bool do_margin;

    assert( rfc_ctx && tp_residue );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    do_margin = rfc_ctx->flags & RFC_FLAGS_ENFORCE_MARGIN;

    if( do_margin && rfc_ctx->tp && !rfc_ctx->tp_locked )
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
                if( !RFC_tp_add( rfc_ctx, &pt_left ) ) return false;

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
bool RFC_feed_finalize( rfc_ctx_s *rfc_ctx )
{
    rfc_value_tuple_s *tp_interim = NULL;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
    
    /* Adjust residue: Incorporate interim turning point */
    if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
    {
        tp_interim = &rfc_ctx->residue[rfc_ctx->residue_cnt];
        rfc_ctx->residue_cnt++;

        rfc_ctx->state = RFC_STATE_BUSY;
    }

#if RFC_TP_SUPPORT
    /* Finalize turning point storage */
    if( !RFC_feed_finalize_tp( rfc_ctx, tp_interim, rfc_ctx->flags ) )
    {
        return false;
    }
#endif /*RFC_TP_SUPPORT*/

    if( tp_interim )
    {
        int flags = rfc_ctx->flags;
#if !RFC_MINIMAL
        /* New turning point, do LC count */
        RFC_cycle_process_lc( rfc_ctx, rfc_ctx->flags & (RFC_FLAGS_COUNT_LC | RFC_FLAGS_ENFORCE_MARGIN) );
        flags &= ~RFC_FLAGS_COUNT_LC;
#endif /*!RFC_MINIMAL*/

        /* Check once more if a new cycle is closed now */
        RFC_cycle_find( rfc_ctx, flags );
    }

#if RFC_HCM_SUPPORT
    if( !RFC_feed_finalize_hcm( rfc_ctx, rfc_ctx->flags ) )
    {
        return false;
    }
#endif /*RFC_HCM_SUPPORT*/

    rfc_ctx->state = RFC_STATE_FINALIZE;

    return true;
}


#if RFC_TP_SUPPORT
/**
 * @brief         Finalize turning point storage.
 *
 * @param         rfc_ctx     The rainflow context
 * @param[in,out] tp_interim  The interim turning point (or NULL)
 * @param[in]     flags       Only flag RFC_FLAGS_ENFORCE_MARGIN encounters
 *
 * @return        true on success
 */
static
bool RFC_feed_finalize_tp( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *tp_interim, int flags )
{
    bool do_margin;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
    
    /* Finalize turning points storage */
    do_margin = rfc_ctx->flags & RFC_FLAGS_ENFORCE_MARGIN;
    if( do_margin && rfc_ctx->tp && !rfc_ctx->tp_locked )
    {
        rfc_value_tuple_s *pt_right = &rfc_ctx->internal.margin[1];

        if( tp_interim )
        {
            if( rfc_ctx->internal.margin_stage > 0 && tp_interim->value == pt_right->value )
            {
                if( !RFC_tp_add( rfc_ctx, pt_right ) ) return false;
                tp_interim->tp_pos = pt_right->tp_pos;
            }
            else
            {
                if( !RFC_tp_add( rfc_ctx, tp_interim ) ) return false;
                if( !RFC_tp_add( rfc_ctx, pt_right ) )   return false;
            }
        }
        else if( pt_right->pos > 0 )
        {
            if( !RFC_tp_add( rfc_ctx, pt_right ) )   return false;
        }
    }
    else if( tp_interim )
    {
        if( !RFC_tp_add( rfc_ctx, tp_interim ) ) return false;
    }

    /* Lock turning points queue */
    RFC_tp_lock( rfc_ctx, true );
    return true;
}
#endif /*RFC_TP_SUPPORT*/


#if RFC_HCM_SUPPORT
/**
 * @brief      Finalize HCM algorithm, copy residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_feed_finalize_hcm( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

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
                return RFC_error_raise( rfc_ctx, RFC_ERROR_MEMORY );
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
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_ignore( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Include interim turning point */
    return RFC_feed_finalize( rfc_ctx );
}


#if !RFC_MINIMAL
/**
 * @brief      Finalize pending counts, discard residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_discard( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Include interim turning point */
    if( !RFC_feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Finalize pending counts, weight unclosed cycles into RP,
 *             RF-matrix and pseudo damage and discard residue.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  weight   The weight for closed cycles
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_weight_cycles( rfc_ctx_s *rfc_ctx, RFC_counts_type weight, int flags )
{
    RFC_counts_type old_inc;
    
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    old_inc = rfc_ctx->curr_inc;

    /* Count every unclosed cycle with the given weight */
    if( rfc_ctx->residue && rfc_ctx->residue_cnt >= 2 )
    {
        size_t             i;
        int                flags = rfc_ctx->flags;
        rfc_value_tuple_s *from  = rfc_ctx->residue;

        rfc_ctx->curr_inc = weight;

        for( i = 0; i + 1 < rfc_ctx->residue_cnt; i++ )
        {
            rfc_value_tuple_s *to   = from + 1;
            rfc_value_tuple_s *next = ( i + 2 < rfc_ctx->residue_cnt ) ? to + 1 : NULL;

            RFC_cycle_process_counts( rfc_ctx, from, to, next, flags );

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
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_clormann_seeger( rfc_ctx_s *rfc_ctx, int flags )
{
    size_t i;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Include interim turning point */
    if( !RFC_feed_finalize( rfc_ctx ) )
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

                RFC_cycle_process_counts( rfc_ctx, from, to, to + 1, flags );

                /* Remove two inner turning points (idx+1 and idx+2) */
                RFC_residue_remove_item( rfc_ctx, i + 1, 2 );
                rfc_ctx->residue_cnt -= 2;
            }
            else
            {
                i++;
            }
        }
    }

    /* Count remaining unclosed cycles half weighted (all but LC counts) */
    return RFC_finalize_res_weight_cycles( rfc_ctx, rfc_ctx->half_inc, flags );
}


/**
 * @brief      Finalize pending counts, DIN method.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_rp_DIN45667( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Include interim turning point */
    if( !RFC_feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* This approach only affects range pair and level crossing countings */
    if( rfc_ctx->flags & RFC_FLAGS_COUNT_RP ) 
    {
        while( rfc_ctx->residue_cnt >= 2 )
        {
            /* Left hand slope to compare */
            rfc_value_tuple_s *from_i       = &rfc_ctx->residue[0];
            rfc_value_tuple_s *to_i         = &rfc_ctx->residue[1];
            int                from_class_i = QUANTIZE( rfc_ctx, from_i->value );
            int                to_class_i   = QUANTIZE( rfc_ctx, to_i->value );
            int                srange_i     = to_class_i - from_class_i;
            size_t             j;

            /* Watch all adjacent slopes */
            for( j = 1; j < rfc_ctx->residue_cnt; j += 2 )
            {
                /* Right hand slopes to compare, all are adjacent to the given "left hand slope" */
                rfc_value_tuple_s *from_j       = &rfc_ctx->residue[j];
                rfc_value_tuple_s *to_j         = &rfc_ctx->residue[j+1];
                int                from_class_j = QUANTIZE( rfc_ctx, from_j->value );
                int                to_class_j   = QUANTIZE( rfc_ctx, to_j->value );
                int                srange_j     = to_class_j - from_class_j;

                /* Matching range found */
                if( srange_i == -srange_j )
                {
                    /* Do the countings for the matching slope */
                    RFC_cycle_process_counts( rfc_ctx, from_j, to_j, to_j + 1, flags );

                    /* Remove "left hand slope" */
                    RFC_residue_remove_item( rfc_ctx, /*index*/ j, /*count*/ 2 );
                }
            }

            /* Remove first point */
            RFC_residue_remove_item( rfc_ctx, /*index*/ 0, /*count*/ 1 );
        }
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Finalize pending counts, repeated residue method.
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    The flags
 *
 * @return     true on success
 */
static
bool RFC_finalize_res_repeated( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Don't include interim turning point! */

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
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
            int                      old_flags = rfc_ctx->flags;
            const rfc_value_tuple_s *from      = rfc_ctx->residue;
                  rfc_value_tuple_s *to        = residue;

            while( n-- )
            {
                *to++ = *from++;
            }

            rfc_ctx->flags = flags;
            /* Feed again with the copy */
            ok = RFC_feed_tuple( rfc_ctx, residue, cnt );
            rfc_ctx->flags = old_flags;

            /* Free temporary residue */
            rfc_ctx->mem_alloc( residue, 0, 0, RFC_MEM_AIM_TEMP );

            if( !ok )
            {
                return false;
            }
        }
        else
        {
            return RFC_error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }
    }

    /* Include interim turning point */
    if( !RFC_feed_finalize( rfc_ctx ) )
    {
        return false;
    }

    /* Empty residue */
    rfc_ctx->residue_cnt = 0;

    return true;
}


/**
 * @brief      Backup/restore of residue
 *
 * @param      rfc_ctx      The rainflow context
 * @param[out] residue      The copy of the current residue
 * @param[out] residue_cap  The capacity of the given residue
 * @param[out] residue_cnt  The number of points in the given residue
 * @param[in]  restore      true->restore, false->backup
 *
 * @return     true on success
 */
static
bool RFC_residue_exchange( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s **residue, size_t *residue_cap, size_t *residue_cnt, bool restore )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
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
            return RFC_error_raise( rfc_ctx, RFC_ERROR_MEMORY );
        }

        memcpy( *residue, rfc_ctx->residue, *residue_cnt * sizeof(rfc_value_tuple_s) );
    }
    else
    {
        /* Restore */

        /* Release residue */
        RFC_mem_alloc( rfc_ctx->residue, /*num*/ 0, /*size*/ 0, RFC_MEM_AIM_TEMP );

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
 * @param[in]  index    The item position in residue, base 0
 * @param[in]  count    The number of elements to remove
 */
static
void RFC_residue_remove_item( rfc_ctx_s *rfc_ctx, size_t index, size_t count )
{
    size_t  from = index + count,
            to   = index, 
            end;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
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
 * @brief      Calculate pseudo damage for one cycle with given amplitude Sa
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  Sa       The amplitude
 *
 * @return     Pseudo damage
 */
static
double RFC_damage_calc_amplitude( rfc_ctx_s *rfc_ctx, double Sa )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    do {
        /* Constants for the Woehler curve */
        const double SD_log = log(rfc_ctx->wl_sd);
        const double ND_log = log(rfc_ctx->wl_nd);
        const double k      = rfc_ctx->wl_k;
#if !RFC_MINIMAL
        const double k2 = rfc_ctx->wl_k2;
#endif /*!RFC_MINIMAL*/

        /* Pseudo damage */
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
                if( Sa > rfc_ctx->wl_sd )
                {
                    D = exp( fabs(k)  * ( log(Sa) - SD_log ) - ND_log );
                }
                else
                {
                    D = exp( fabs(k2) * ( log(Sa) - SD_log ) - ND_log );
                }
            }
#else /*RFC_MINIMAL*/
            D = exp( fabs(k)  * ( log(Sa) - SD_log ) - ND_log );
#endif /*!RFC_MINIMAL*/
        }
        else
        {
            assert( false );
        }

        return D;
        
    } while(0);
}


#if RFC_AT_SUPPORT
/**
 * @brief      Calculate the normalized mean load (Sa=1) for a given load ratio
 *             "R"
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  R        Load ratio (Su/So)
 *
 * @return     Normalized mean load "Sm_norm"
 */
static
double RFC_at_R_to_Sm_norm( rfc_ctx_s *rfc_ctx, double R )
{
    double Sm_norm;

    if( isinf( R ) )
    {
        Sm_norm = -1.0;
    }
    else
    {
        /* Su      = Sm - Sa */
        /* So      = Sm + Sa */
        /* R       = Su / So */
        /* Sm      = Sa * ( 1 + R ) / ( 1 - R ) */
        /* Sm_norm = Sm / Sa */

        Sm_norm = ( 1 + R ) / ( 1 - R );
    }

    return Sm_norm;
}


/**
 * @brief      Calculate the influence of (normalized, Sa=1) mean load on
 *             fatigue strength (Haigh-diagram)
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  Sm_norm  Mean load, normalized (Sm/Sa)
 *
 * @return     Alleviation factor on reference curve
 */
static
double RFC_at_alleviation( rfc_ctx_s *rfc_ctx, double Sm_norm )
{
    assert( rfc_ctx );

    if( !rfc_ctx->at.count )
    {
        /* No reference curve given, no transformation */
        return 1.0;
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
            return Sa_[0];  /* Clip to first point */
        }
        else if( Sm_norm >= Sm_[count-1] / Sa_[count-1] )
        {
            /* Last segment */
            return Sa_[count-1];  /* Clip to last point */
        }
        else
        {
            unsigned i;

            /* Select correct segment */
            for( i = 1; i < count; i++ )
            {
                assert( Sa_[i-1] > 0.0 && Sa_[i] > 0.0 && Sm_[i-1] <= Sm_[i] );

                if( Sm_norm > Sm_[i-1] / Sa_[i-1] && Sm_norm <= Sm_[i] / Sa_[i] )
                {
                    /* Intersection of R slope and M slope */
                    double M_signed    = ( Sa_[i] - Sa_[i-1] ) / ( Sm_[i] - Sm_[i-1] );
                    double alleviation = ( Sa_[i-1] - M_signed * Sm_[i-1] ) / ( 1.0 - M_signed * Sm_norm );

                    return alleviation;
                }
            }
        }

        assert( false );
        return 1.0;
    }
}
#endif /*RFC_AT_SUPPORT*/


/**
 * @brief      Calculate fictive damage for one closed (full) cycle.
 *
 * @param      rfc_ctx     The rainflow context
 * @param[in]  class_from  The starting class
 * @param[in]  class_to    The ending class
 *
 * @return     Pseudo damage value for the closed cycle
 */
static
double RFC_damage_calc( rfc_ctx_s *rfc_ctx, unsigned class_from, unsigned class_to )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

#if RFC_DAMAGE_FAST
    if( rfc_ctx->damage_lut )
    {
        return RFC_damage_calc_fast( rfc_ctx, class_from, class_to );
    }
#endif /*RFC_DAMAGE_FAST*/

#if RFC_USE_DELEGATES
    if( rfc_ctx->damage_calc_fcn )
    {
        return rfc_ctx->damage_calc_fcn( rfc_ctx, class_from, class_to );
    }
#endif /*RFC_USE_DELEGATES*/


    if( class_from != class_to )
    {
#if RFC_MINIMAL
        double range     = (double)rfc_ctx->class_width * abs( (int)class_to - (int)class_from );
        double amplitude = range / 2.0;

        return RFC_damage_calc_amplitude( rfc_ctx, amplitude );
#else /*!RFC_MINIMAL*/
        double Sa_i   = fabs( (int)class_from - (int)class_to ) / 2.0 * rfc_ctx->class_width;
        double Sm_i   =     ( (int)class_from + (int)class_to ) / 2.0 * rfc_ctx->class_width + rfc_ctx->class_offset;

        if( Sa_i > 0.0 )
        {
#if RFC_AT_SUPPORT
            /* Calculate transformation factor with normalized mean value */
            Sa_i = RFC_at_transform( rfc_ctx, Sa_i, Sm_i );
#endif /*RFC_AT_SUPPORT*/

            return RFC_damage_calc_amplitude( rfc_ctx, Sa_i );
        }
#endif /*RFC_MINIMAL*/
    }

    return 0.0;
}


#if RFC_DAMAGE_FAST
/**
 * @brief      Initialize a look-up table of damages for closed cycles. In this
 *             implementation the midrange doesn't matter!
 *
 * @param      rfc_ctx  The rainflow context
 */
static 
void RFC_damage_lut_init( rfc_ctx_s *rfc_ctx )
{
    double   *lut;
    unsigned  from, 
              to;

    assert( rfc_ctx );
    assert( rfc_ctx->state == RFC_STATE_INIT );

    if( rfc_ctx->damage_lut )
    {
        lut = rfc_ctx->damage_lut;
        rfc_ctx->damage_lut = NULL;

        for( from = 0; from < rfc_ctx->class_count; from++ )
        {
            for( to = 0; to < rfc_ctx->class_count; to++ )
            {
                lut[from * rfc_ctx->class_count + to] = RFC_damage_calc( rfc_ctx, from, to );
            }
        }

        rfc_ctx->damage_lut = lut;
    }
}


/**
 * @brief      Calculate fictive damage for one closed (full) cycle, using
 *             look-up table.
 *
 * @param      rfc_ctx     The rainflow context
 * @param[in]  class_from  The starting class
 * @param[in]  class_to    The ending class
 *
 * @return     Pseudo damage value for the closed cycle
 */
static
double RFC_damage_calc_fast( rfc_ctx_s *rfc_ctx, unsigned class_from, unsigned class_to )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );
    assert( class_from < rfc_ctx->class_count );
    assert( class_to   < rfc_ctx->class_count );

    if( rfc_ctx->damage_lut )
    {
        return( rfc_ctx->damage_lut[class_from * rfc_ctx->class_count + class_to] );
    }

    return 0.0;
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
rfc_value_tuple_s * RFC_feed_filter_pt( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *pt )
{
    int                 slope;
    RFC_value_type      delta;
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
 * @param[in]  flags    The flags
 */
static
void RFC_cycle_find( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

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
                RFC_cycle_find_4ptm( rfc_ctx, flags );
                break;
#if RFC_HCM_SUPPORT
            case RFC_COUNTING_METHOD_HCM:
                RFC_cycle_find_hcm( rfc_ctx, flags );
                break;
#endif /*RFC_HCM_SUPPORT*/
            case RFC_COUNTING_METHOD_DELEGATED:
                /* FALLTHROUGH */
            default:
                assert( false );
                break;
        }
    }

    if( rfc_ctx->counting_method == RFC_COUNTING_METHOD_NONE || !rfc_ctx->class_count )
    {
        /* Prune residue */
        if( rfc_ctx->residue_cnt > 1 )
        {
            RFC_residue_remove_item( rfc_ctx, /*index*/ 0, rfc_ctx->residue_cnt - 1 );
        }
    }

}
#endif /*!RFC_MINIMAL*/


/**
 * @brief      Rainflow counting core (4-point-method).
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    The flags
 */
static
void RFC_cycle_find_4ptm( rfc_ctx_s *rfc_ctx, int flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    while( rfc_ctx->residue_cnt >= 4 )
    {
        size_t idx = rfc_ctx->residue_cnt - 4;

        RFC_value_type A = rfc_ctx->residue[idx+0].value;
        RFC_value_type B = rfc_ctx->residue[idx+1].value;
        RFC_value_type C = rfc_ctx->residue[idx+2].value;
        RFC_value_type D = rfc_ctx->residue[idx+3].value;

        if( B > C )
        {
            RFC_value_type temp = B;
            B = C;
            C = temp;
        }

        if( A > D )
        {
            RFC_value_type temp = A;
            A = D;
            D = temp;
        }

        if( A <= B && C <= D )
        {
            rfc_value_tuple_s *from = &rfc_ctx->residue[idx+1];
            rfc_value_tuple_s *to   = &rfc_ctx->residue[idx+2];

            /* Closed cycle found, process countings */
            RFC_cycle_process_counts( rfc_ctx, from, to, to + 1, flags );

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
 * @param[in]  flags    The flags
 */
static
void RFC_cycle_find_hcm( rfc_ctx_s *rfc_ctx, int flags )
{
    int IZ, IR;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

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

            if( (K->value - J->value) * (J->value - I->value) >= 0 )
            {
                /* Is no turning point */
                /* This should only may happen, when RFC_FLAGS_ENFORCE_MARGIN is set, 
                   since all values from residue are turning points */
                assert( rfc_ctx->flags & RFC_FLAGS_ENFORCE_MARGIN );
                IZ--;
                /* Test further closed cycles */
                goto label_2;
            }
            else
            {
                /* Is a turning point */
                if( fabs( (double)K->value - (double)J->value ) >= fabs( (double)J->value - (double)I->value) )
                {
                    /* Cycle range is greater or equal to previous, register closed cycle */
                    RFC_cycle_process_counts( rfc_ctx, I, J, NULL, flags );
                    IZ -= 2;
                    /* Test further closed cycles */
                    goto label_2;
                }
            }
        }
        else if( IZ == IR )
        {
            J = &rfc_ctx->internal.hcm.stack[IZ];

            if( ( (double)K->value - (double)J->value ) * (double)J->value >= 0.0 )
            {
                /* Is no turning point */
                IZ--;
                /* Test further closed cycles */
                goto label_2;
            }
            else if( fabs( (double)K->value ) > fabs( (double)J->value ) )
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
        RFC_residue_remove_item( rfc_ctx, /*index*/ 0, 1 );
    }

    /* hcm.IZ and hcm.IR are base 1! */
    rfc_ctx->internal.hcm.IZ = IZ + 1;
    rfc_ctx->internal.hcm.IR = IR + 1;
}
#endif /*RFC_HCM_SUPPORT*/


#if !RFC_MINIMAL
/**
 * @brief      Processes LC count (level crossing)
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  flags    Control flags
 */
static
void RFC_cycle_process_lc( rfc_ctx_s *rfc_ctx, int flags )
{
    size_t n = rfc_ctx->residue_cnt;

    if( n > 1 && (flags & RFC_FLAGS_COUNT_LC) )
    {
        /* Do the level crossing counting */
        RFC_cycle_process_counts( rfc_ctx, &rfc_ctx->residue[n-2], &rfc_ctx->residue[n-1], NULL, flags & (RFC_FLAGS_COUNT_LC | RFC_FLAGS_ENFORCE_MARGIN) );
    }
}
#endif /*!RFC_MINIMAL*/


/**
 * @brief         Processes counts on a closing cycle
 *
 * @param         rfc_ctx  The rainflow context
 * @param[in,out] from     The starting data point
 * @param[in,out] to       The ending data point
 * @param[in,out] next     The point next after "to"
 * @param[in]     flags    Control flags
 */
static
void RFC_cycle_process_counts( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, int flags )
{
    unsigned class_from, class_to;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
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

    /* Quantize "from" */
    class_from = QUANTIZE( rfc_ctx, from->value );

    if( class_from >= rfc_ctx->class_count ) class_from = rfc_ctx->class_count - 1;

    /* Quantize "to" */
    class_to = QUANTIZE( rfc_ctx, to->value );

    if( class_to >= rfc_ctx->class_count ) class_to = rfc_ctx->class_count - 1;
    
    /* class_from and class_to are base 0 now */

    /* Do several counts, according to "flags" */
    if( class_from != class_to )
    {
        /* Cumulate pseudo damage */
        if( flags & RFC_FLAGS_COUNT_DAMAGE )
        {
            double damage = RFC_damage_calc( rfc_ctx, class_from, class_to );

            /* Adding damage due to current cycle weight */
            rfc_ctx->pseudo_damage += damage * rfc_ctx->curr_inc / rfc_ctx->full_inc;
        }

        /* Rainflow matrix */
        if( rfc_ctx->matrix && ( flags & RFC_FLAGS_COUNT_MATRIX ) )
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
            
            assert( rfc_ctx->matrix[idx] <= RFC_COUNTS_LIMIT );
            rfc_ctx->matrix[idx] += rfc_ctx->curr_inc;
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
            
            if( class_from < class_to && ( flags & RFC_FLAGS_COUNT_LC_UP ) )
            {
                /* Count rising slopes */
                unsigned idx;
                for( idx = class_from; idx < class_to; idx++ )
                {
                    assert( rfc_ctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfc_ctx->lc[idx] += rfc_ctx->curr_inc;
                }
            }
            else if( class_to < class_from && ( flags & RFC_FLAGS_COUNT_LC_DN ) )
            {
                /* Count falling slopes */
                unsigned idx;
                for( idx = class_to; idx < class_from; idx++ )
                {
                    assert( rfc_ctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfc_ctx->lc[idx] += rfc_ctx->curr_inc;
                }
            }
        }

#if RFC_DH_SUPPORT
        if( flags & RFC_FLAGS_COUNT_DH )
        {
            /* "Spread" damage over turning points (tp) and damage history (dh) */
            RFC_spread_damage( rfc_ctx, from, to, next, flags );
        }
#endif /*RFC_DH_SUPPORT*/

#endif /*!RFC_MINIMAL*/
    }
}


#if RFC_TP_SUPPORT
/**
 * @brief         Append one turning point to the queue
 *
 * @param         rfc_ctx  The rainflow context
 * @param[in,out] tp       New turning points
 *
 * @return        true on success
 */
static
bool RFC_tp_add( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *tp )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Add new turning point */

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->tp_add_fcn )
    {
        /* Add turning point */
        return rfc_ctx->tp_add_fcn( rfc_ctx, tp );
    }
    else
#endif /*RFC_USE_DELEGATES*/
    {
        if( tp && rfc_ctx->tp && !rfc_ctx->tp_locked )
        {
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
                    return RFC_error_raise( rfc_ctx, RFC_ERROR_MEMORY );
                }
            }
            /* Append turning point */
            rfc_ctx->tp[ rfc_ctx->tp_cnt++ ] = *tp;
            tp->tp_pos = rfc_ctx->tp_cnt; /* Turning point get its position index in tp, base 1 */
        }

        if( rfc_ctx->flags & RFC_FLAGS_TPAUTOPRUNE && rfc_ctx->tp_cnt > rfc_ctx->tp_prune_threshold )
        {
            return RFC_tp_prune( rfc_ctx, rfc_ctx->tp_prune_size, RFC_FLAGS_TPPRUNE_PRESERVE_POS );
        }

        return true;
    }
}


/**
 * @brief      Lock turning points queue
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  do_lock  Turning point storage will be locked, if true
 */
static 
void RFC_tp_lock( rfc_ctx_s *rfc_ctx, bool do_lock )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    rfc_ctx->tp_locked = do_lock;
}


/**
 * @brief      Refeed all values from turning points history
 *
 * @param      rfc_ctx          The rainflow context
 * @param[in]  new_hysteresis   The new hysteresis
 * @param[in]  new_class_param  The new class parameters
 */
static
void RFC_tp_refeed( rfc_ctx_s *rfc_ctx, RFC_value_type new_hysteresis, const rfc_class_param_s *new_class_param )
{
    rfc_value_tuple_s *ctx_tp;
    size_t             ctx_tp_cnt,
                       i;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Get the current turning points */
    ctx_tp     = rfc_ctx->tp;
    ctx_tp_cnt = rfc_ctx->tp_cnt;

    /* RFC_reset sets rfc_ctx->tp_cnt to zero, but turning points are still available */
    RFC_reset( rfc_ctx );

    /* Class parameters may change, new hysteresis must be greater! */
    if( new_class_param )
    {
        assert( new_hysteresis >= rfc_ctx->hysteresis );
        rfc_ctx->hysteresis = new_hysteresis;
        RFC_class_param_set( rfc_ctx, new_class_param );
#if RFC_DAMAGE_FAST
        RFC_damage_lut_init( rfc_ctx );
#endif /*RFC_DH_SUPPORT*/
    }

    for( i = 0; i < ctx_tp_cnt; i++, ctx_tp++ )
    {
        ctx_tp->class = QUANTIZE( rfc_ctx, ctx_tp->value );
    }

    RFC_feed_tuple( rfc_ctx, ctx_tp, ctx_tp_cnt );

}
#endif /*RFC_TP_SUPPORT*/


#if RFC_DH_SUPPORT
/**
 * @brief         Spread damage over damage history
 *
 * @param         rfc_ctx  The rainflow context
 * @param[in,out] from     The start turning point
 * @param[in,out] to       The end turning point
 * @param[in,out] next     The next turning point after @a to
 * @param[in]     flags    The flags
 */
static 
void RFC_spread_damage( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s *from, 
                                            rfc_value_tuple_s *to, 
                                            rfc_value_tuple_s *next, int flags )
{
    double damage = 0.0;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINALIZE );
    assert( from && to );

    switch( rfc_ctx->spread_damage_method )
    {
        case RFC_SD_NONE:
            break;
        case RFC_SD_HALF_23:
            damage = RFC_damage_calc( rfc_ctx, from->class, to->class );
#if RFC_TP_SUPPORT && RFC_DH_SUPPORT
            from->damage += damage / 2.0;
            to->damage   += damage / 2.0;
#endif /*RFC_TP_SUPPORT && RFC_DH_SUPPORT*/
            break;
        case RFC_SD_FULL_P2:
            damage = RFC_damage_calc( rfc_ctx, from->class, to->class );
#if RFC_TP_SUPPORT && RFC_DH_SUPPORT
            from->damage += damage;
#endif /*RFC_TP_SUPPORT && RFC_DH_SUPPORT*/
            break;
        case RFC_SD_FULL_P3:
            damage = RFC_damage_calc( rfc_ctx, from->class, to->class );
#if RFC_TP_SUPPORT && RFC_DH_SUPPORT
            to->damage += damage;
#endif /*RFC_TP_SUPPORT && RFC_DH_SUPPORT*/
            break;
        case RFC_SD_RAMP_AMPLITUDE_23:
        case RFC_SD_RAMP_DAMAGE_23:
        case RFC_SD_RAMP_AMPLITUDE_24:
        case RFC_SD_RAMP_DAMAGE_24:
        {
#if RFC_TP_SUPPORT
            size_t  i,
                    start, end, width, 
                    tp_start, tp_end;
            int     range;

            /* Care about possible wrapping, caused by RFC_RES_REPEATED:
                0         1
                01234567890123 (14 points)
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

            /* Absolute position (input stream) */
            start    = from->pos;
            end      = to->pos;
            end     += ( start >= end ) ? rfc_ctx->internal.pos : 0;
            width    = end - start;
            /* Position in turning point storage */
            tp_start = from->tp_pos;
            tp_end   = to->tp_pos;
            tp_end  += ( tp_start >= tp_end ) ? rfc_ctx->tp_cnt : 0;

            range    = abs( (int)to->class - (int)from->class );
            damage   = 0.0;

            /* Iterate over turning points */
            for( i = tp_start; i <= tp_end; i++ )
            {
                size_t tp_pos, pos;
                double weight, new_damage;

                tp_pos = i % rfc_ctx->tp_cnt;
                pos    = rfc_ctx->tp[tp_pos].pos;
                pos   += ( start >= pos ) ? rfc_ctx->internal.pos : 0;
                weight = (double)( pos - start ) / width;

                switch( rfc_ctx->spread_damage_method )
                {
                    case RFC_SD_RAMP_AMPLITUDE_23:
                    case RFC_SD_RAMP_AMPLITUDE_24:
                        damage = RFC_damage_calc( rfc_ctx, from->class, to->class );
                        break;
                    case RFC_SD_RAMP_DAMAGE_23:
                    case RFC_SD_RAMP_DAMAGE_24:
                        damage = RFC_damage_calc( rfc_ctx, from->class, to->class );
                        break;
                }

                if( new_damage > damage )
                {
                    rfc_ctx->tp[tp_pos].damage += new_damage - damage;
                    damage = new_damage;
                }
            }
#endif /*RFC_TP_SUPPORT*/
        }
        break;

        case RFC_SD_TRANSIENT_23:
            break;
        case RFC_SD_TRANSIENT_23c:
            break;

        default:
            assert( false );
            break;
    }
}
#endif /*RFC_TP_SUPPORT*/


/**
 * @brief      Returns the unsigned difference of two values, sign optionally
 *             returned as -1 or 1.
 *
 * @param[in]  from      Left hand value
 * @param[in]  to        Right hand value
 * @param[out] sign_ptr  Pointer to catch sign (may be NULL)
 *
 * @return     Returns the absolute difference of given values
 */
static
RFC_value_type value_delta( RFC_value_type from, RFC_value_type to, int *sign_ptr )
{
    double delta = (double)to - (double)from;

    if( sign_ptr )
    {
        *sign_ptr = ( delta < 0.0 ) ? -1 : 1;
    }

    return (RFC_value_type)fabs( delta );
}


/**
 * @brief      Raises an error
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  error    The error identifier
 *
 * @return     Always false
 */
static
bool RFC_error_raise( rfc_ctx_s *rfc_ctx, int error )
{
    assert( rfc_ctx );

    rfc_ctx->state = RFC_STATE_ERROR;
    rfc_ctx->error = error;

    return false;
}


/**
 * @brief      (Re-)Allocate or free memory
 *
 * @param      ptr   Previous data pointer, or NULL, if unset
 * @param[in]  num   The number of elements
 * @param[in]  size  The size of one element in bytes
 * @param[in]  aim   The aim
 *
 * @return     New memory pointer or NULL if either num or size is 0
 */
static
void * RFC_mem_alloc( void *ptr, size_t num, size_t size, int aim )
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

#if !RFC_MINIMAL
/**
 * @brief      Set class parameter.
 *
 * @param      rfc_ctx      The rainflow context
 * @param[in]  class_param  The new class parameter
 */
static
void RFC_class_param_set( rfc_ctx_s *rfc_ctx, const rfc_class_param_s *class_param )
{
    assert( rfc_ctx && class_param );
    assert( rfc_ctx->state == RFC_STATE_INIT );

    assert( class_param->count > 0 );
    assert( class_param->width > 0.0 );

    rfc_ctx->class_count  = class_param->count;
    rfc_ctx->class_width  = class_param->width;
    rfc_ctx->class_offset = class_param->offset;
}


/**
 * @brief      Get class parameter.
 *
 * @param      rfc_ctx      The rainflow context
 * @param[in]  class_param  The class parameter
 */
static
void RFC_class_param_get( rfc_ctx_s *rfc_ctx, rfc_class_param_s *class_param )
{
    assert( rfc_ctx && class_param );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    class_param->count  = rfc_ctx->class_count;
    class_param->width  = rfc_ctx->class_width;
    class_param->offset = rfc_ctx->class_offset;
}
#endif /*!RFC_MINIMAL*/

/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/

#if !RFC_MINIMAL
bool RFC_lc_from_matrix     ( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags );
bool RFC_lc_from_residue    ( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags );
bool RFC_rp_from_matrix     ( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags );

/**
 * Calculate level crossing counts from rainflow matrix, write results to
 * buffer.
 *
 * @param      rfc_ctx     The rainflow context
 * @param      buffer      The buffer for LC histogram
 * @param[in]  buffer_cap  The buffer capability
 * @param[in]  flags       The flags
 *
 * @return     true on success
 */
bool RFC_lc_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags )
{
    unsigned i, j, k;
    bool     up     = flags & RFC_FLAGS_COUNT_LC_UP;
    bool     dn     = flags & RFC_FLAGS_COUNT_LC_DN;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    if( !buffer || buffer_cap < rfc_ctx->class_count ) return false;

    for( i = 0; i < rfc_ctx->class_count; i++ ) 
    {
        RFC_counts_type counts = (RFC_counts_type)0;

        for( j = i; j < rfc_ctx->class_count; j++ )  /* To */
        {
            for( k = 0; k < i; k++ )                 /* From */
            {
                /* Count rising slopes */
                assert( counts < RFC_COUNTS_LIMIT - rfc_ctx->matrix[ k * rfc_ctx->class_count + j ] );
                if( up )
                {
                    counts += rfc_ctx->matrix[ k * rfc_ctx->class_count + j ];
                }

                /* Count falling slopes */
                assert( counts < RFC_COUNTS_LIMIT - rfc_ctx->matrix[ j * rfc_ctx->class_count + k ] );
                if( dn )
                {
                    counts += rfc_ctx->matrix[ j * rfc_ctx->class_count + k ];
                }
            }
        }

        buffer[i] = counts;
    }

    return true;
}


/**
 * Calculate level crossing counts from rainflow matrix, write results to
 * buffer.
 *
 * @param      rfc_ctx     The rainflow context
 * @param[out] buffer      The buffer for LC histogram
 * @param[in]  buffer_cap  The buffer capability
 * @param[in]  flags       The flags
 *
 * @return     true on success
 */

bool RFC_lc_from_residue( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags )
{
    unsigned i;
    bool     up     = flags & RFC_FLAGS_COUNT_LC_UP;
    bool     dn     = flags & RFC_FLAGS_COUNT_LC_DN;
    const rfc_value_tuple_s *from = rfc_ctx->residue;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    if( !buffer || buffer_cap < rfc_ctx->class_count ) return false;

    for( i = 1; i < rfc_ctx->residue_cnt; i++ ) 
    {
        const rfc_value_tuple_s *to         = from + 1;
              unsigned           class_from = from->class;
              unsigned           class_to   = to->class;

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
                assert( buffer[idx] <= RFC_COUNTS_LIMIT );
                buffer[idx] += rfc_ctx->full_inc;
            }
        }
        else if( class_to < class_from && dn )
        {
            /* Count falling slopes */
            unsigned idx;
            for( idx = class_to; idx < class_from; idx++ )
            {
                assert( buffer[idx] <= RFC_COUNTS_LIMIT );
                buffer[idx] += rfc_ctx->full_inc;
            }
        }

        from = to;
    }

    return true;
}


/**
 * Calculate range pair counts from rainflow matrix, write results to buffer.
 *
 * @param      rfc_ctx     The rainflow context
 * @param[out] buffer      The buffer for RP histogram
 * @param[in]  buffer_cap  The buffer capability
 * @param[in]  flags       The flags
 *
 * @return     true on success
 */
bool RFC_rp_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_cap, int flags )
{
    unsigned i, j;
    size_t   maxcnt = buffer_cap / sizeof(RFC_counts_type);

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );

    if( !buffer || buffer_cap < rfc_ctx->class_count ) return false;

    buffer[0] = 0;

    for( i = 1; i < rfc_ctx->class_count; i++ ) 
    {
        RFC_counts_type counts = (RFC_counts_type)0;

        for( j = i; j < rfc_ctx->class_count; j++ )
        {
            /* Count rising and falling slopes */
            assert( counts < RFC_COUNTS_LIMIT - rfc_ctx->matrix[ i * rfc_ctx->class_count + j ] 
                                              - rfc_ctx->matrix[ j * rfc_ctx->class_count + i ] );
            counts += rfc_ctx->matrix[ i * rfc_ctx->class_count + j ];
            counts += rfc_ctx->matrix[ j * rfc_ctx->class_count + i ];
        }

        buffer[i > maxcnt ? maxcnt : i] = counts;
    }

    return true;
}

#endif /*!RFC_MINIMAL*/


/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


#if MATLAB_MEX_FILE

#if 0
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    mexErrMsgTxt( "Unsupported configuration!" );
}
#else

/**
 * MATLAB wrapper for the rainflow algorithm
 */
static
void mexRainflow( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
#if !RFC_MINIMAL
    if( nrhs != 8 )
    {
        mexErrMsgTxt( "Function needs exact 7 arguments!" );
#else /*RFC_MINIMAL*/
    if( nrhs != 6 )
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
        const mxArray  *mxResidualMethod = prhs[5];
#if !RFC_MINIMAL
        const mxArray  *mxEnforceMargin  = prhs[6];
        const mxArray  *mxUseHCM         = prhs[7];
#endif /*!RFC_MINIMAL*/        

        RFC_value_type *buffer           = NULL;
        double         *data             = mxGetPr( mxData );
        size_t          data_len         = mxGetNumberOfElements( mxData );
        unsigned        class_count      = (unsigned)( mxGetScalar( mxClassCount ) + 0.5 );
        double          class_width      = mxGetScalar( mxClassWidth );
        double          class_offset     = mxGetScalar( mxClassOffset );
        double          hysteresis       = mxGetScalar( mxHysteresis );
        int             residual_method  = (int)( mxGetScalar( mxResidualMethod ) + 0.5 );
#if !RFC_MINIMAL
        int             enforce_margin   = (int)mxGetScalar( mxEnforceMargin );
        int             use_hcm          = (int)mxGetScalar( mxUseHCM );
#endif /*!RFC_MINIMAL*/
        size_t          i;
        bool            ok;

        mxAssert( residual_method >= RFC_RES_NONE && residual_method <= RFC_RES_RP_DIN45667, 
                  "Invalid residual method!" );

        ok = RFC_init( &rfc_ctx, 
                       class_count, (RFC_value_type)class_width, (RFC_value_type)class_offset, 
                       (RFC_value_type)hysteresis );
#if RFC_TP_SUPPORT
        if( ok )
        {
            rfc_ctx.tp_cap = 128;
            rfc_ctx.tp     = (rfc_value_tuple_s*)RFC_mem_alloc( NULL, rfc_ctx.tp_cap, 
                                                                sizeof(rfc_value_tuple_s), RFC_MEM_AIM_TP );

            ok = RFC_tp_init( &rfc_ctx, rfc_ctx.tp, rfc_ctx.tp_cap, /* is_static */ true );
        }
#endif /*RFC_TP_SUPPORT*/

        if( !ok )
        {
            mexErrMsgTxt( "Error during initialization!" );
        }

        /* Cast values from double type to RFC_value_type */ 
        if( sizeof( RFC_value_type ) != sizeof(double) && data_len )  /* maybe unsafe! */
        {
            buffer = (RFC_value_type *)RFC_mem_alloc( NULL, data_len, 
                                                      sizeof(RFC_value_type), RFC_MEM_AIM_TEMP );

            if( !buffer )
            {
                RFC_deinit( &rfc_ctx );
                mexErrMsgTxt( "Error during initialization!" );
            }

            for( i = 0; i < data_len; i++ )
            {
                buffer[i] = (RFC_value_type)data[i];
            }
        }
        else buffer = (RFC_value_type*)data;

        /* Rainflow counting */

#if !RFC_MINIMAL
        /* Setup */
        rfc_ctx.flags           |= enforce_margin ? RFC_FLAGS_ENFORCE_MARGIN : 0;
#endif /*!RFC_MINIMAL*/
#if RFC_HCM_SUPPORT
        rfc_ctx.counting_method  = use_hcm ? RFC_COUNTING_METHOD_HCM : RFC_COUNTING_METHOD_4PTM;
#else /*!RFC_HCM_SUPPORT*/
#if !RFC_MINIMAL
        rfc_ctx.counting_method  = RFC_COUNTING_METHOD_4PTM;
#endif /*!RFC_MINIMAL*/
#endif /*RFC_HCM_SUPPORT*/
        RFC_feed( &rfc_ctx, buffer, data_len );
        RFC_finalize( &rfc_ctx, residual_method );

        /* Free temporary buffer (cast) */
        if( (void*)buffer != (void*)data )
        {
            RFC_mem_alloc( buffer, 0, 0, RFC_MEM_AIM_TEMP );
            buffer = NULL;
        }

        /* Return results */
        if( plhs )
        {
            /* Pseudo damage */
            plhs[0] = mxCreateDoubleScalar( rfc_ctx.pseudo_damage );

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
            if( nlhs > 2 && rfc_ctx.matrix )
            {
                mxArray* matrix = mxCreateDoubleMatrix( class_count, class_count, mxREAL );
                if( matrix )
                {
                    double *ptr = mxGetPr(matrix);
                    size_t from, to;
                    for( to = 0; to < class_count; to++ )
                    {
                        for( from = 0; from < class_count; from++ )
                        {
                            *ptr++ = (double)rfc_ctx.matrix[ from * class_count + to ] / rfc_ctx.full_inc;
                        }
                    }
                    plhs[2] = matrix;
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
                mxArray* tp = mxCreateDoubleMatrix( rfc_ctx.tp_cnt, 2, mxREAL );
                if( tp )
                {
                    size_t i;
                    double *idx = mxGetPr(tp) + 0;
                    double *val = mxGetPr(tp) + rfc_ctx.tp_cnt;

                    for( i = 0; i < rfc_ctx.tp_cnt; i++ )
                    {
                        *val++ = (double)rfc_ctx.tp[i].value;
                        *idx++ = (double)rfc_ctx.tp[i].pos;
                    }
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
#endif /*0~1*/


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

        mxAssert( RFC_init( &ctx, 0 /*class_count*/, 0.0 /*class_width*/, 0.0 /*class_offset*/, 0.0 /*hysteresis*/ ),
                  "RFC initialization error!" );
        mxAssert( RFC_at_init( &ctx, NULL /*Sa*/, NULL /*Sm*/, 0 /*count*/, M, target /*Sm_rig*/, target /*R_rig*/, target_is_R, false /*symmetric*/ ), 
                  "RFC initialization error!" );

        for( n = 0; n < cnt; n++ )
        {
            double Sa_n = mxGetPr( mxSa )[n];
            double Sm_n = mxGetPr( mxSm )[n];

            mxGetPr( mxResult )[n] = RFC_at_transform( &ctx, Sa_n, Sm_n );
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

        tp = calloc( mxGetNumberOfElements( mxData ), sizeof(tp[0]) );
        if( !tp )
        {
            mexErrMsgTxt( "Memory allocation error!" );
        }

        if( !RFC_init( &ctx, 0 /*class_count*/, 0.0 /*class_width*/, 0.0 /*class_offset*/, hysteresis ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC init!" );
        }

        if( !RFC_tp_init( &ctx, tp, mxGetNumberOfElements( mxData ), true /*tp_is_static*/ ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC tp init!" );
        }

        RFC_flags_set( &ctx, enforce_margin ? RFC_FLAGS_ENFORCE_MARGIN : 0 );

        if( !RFC_feed( &ctx, mxGetPr( mxData ), mxGetNumberOfElements( mxData ) ) )
        {
            free(tp);
            mexErrMsgTxt( "Error on RFC feed!" );
        }

        if( !RFC_finalize_res_ignore( &ctx, ctx.flags ) )
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
#endif /*!RFC_MINIMAL*/
    
#endif /*MATLAB_MEX_FILE*/
