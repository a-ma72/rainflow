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
 * []  "Zaehlverfahren und Lastannahme in der Betriebsfestigkeit";
 *     Michael Koehler, Sven Jenne / Kurt Poetter, Harald Zenner; Springer-Verlag Berlin Heidelberg 2012
 *
 *
 *================================================================================
 * BSD 2-Clause License
 * 
 * Copyright (c) 2023, Andras Martin
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

static char* __rfc_core_version__ = RFC_CORE_VERSION;

#ifndef CALLOC
#define CALLOC calloc
#endif
#ifndef REALLOC
#define REALLOC realloc
#endif
#ifndef FREE
#define FREE free
#endif


/* Core functions */
#define cycle_find          cycle_find_4ptm
static bool                 feed_once                       (       rfc_ctx_s *, const rfc_value_tuple_s* tp, rfc_flags_e flags );
static bool                 feed_finalize                   (       rfc_ctx_s * );
static rfc_value_tuple_s *  feed_filter_pt                  (       rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 cycle_find_4ptm                 (       rfc_ctx_s *, rfc_flags_e flags );
static void                 cycle_process_counts            (       rfc_ctx_s *, rfc_value_tuple_s *from, rfc_value_tuple_s *to, rfc_value_tuple_s *next, rfc_flags_e flags );
/* Methods on residue */
static bool                 finalize_res_ignore             (       rfc_ctx_s *, rfc_flags_e flags );
static bool                 finalize_res_no_finalize        (       rfc_ctx_s *, rfc_flags_e flags );
static void                 residue_remove_item             (       rfc_ctx_s *, size_t index, size_t count );
/* Memory allocator */
static void *               mem_alloc                       ( void *ptr, size_t num, size_t size, int aim );
/* Other */
static bool                 damage_calc_amplitude           (       rfc_ctx_s *, double Sa, double *damage );
static bool                 damage_calc                     (       rfc_ctx_s *, unsigned class_from, unsigned class_to, double *damage, double *Sa_ret );
static bool                 error_raise                     (       rfc_ctx_s *, rfc_error_e );
static rfc_value_t          value_delta                     (       rfc_ctx_s *, const rfc_value_tuple_s* pt_from, const rfc_value_tuple_s* pt_to, int *sign_ptr );


#define QUANTIZE( r, v )    ( (r)->class_count ? (unsigned)( ((v) - (r)->class_offset) / (r)->class_width ) : 0 )
#define AMPLITUDE( r, i )   ( (r)->class_count ? ( (double)(r)->class_width * (i) / 2 ) : 0.0 )
#define CLASS_MEAN( r, c )  ( (r)->class_count ? ( (double)(r)->class_width * (0.5 + (c)) + (r)->class_offset ) : 0.0 )
#define CLASS_UPPER( r, c ) ( (r)->class_count ? ( (double)(r)->class_width * (1.0 + (c)) + (r)->class_offset ) : 0.0 )
#define NUMEL( x )          ( sizeof(x) / sizeof(*(x)) )
#define MAT_OFFS( i, j )    ( (i) * class_count + (j) )

#define RFC_CTX_CHECK_AND_ASSIGN                                                    \
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;                                           \
                                                                                    \
    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )                         \
    {                                                                               \
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );                            \
    }                                                                               \



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
        flags                               = RFC_FLAGS_COUNT_RFM            | 
                                              RFC_FLAGS_COUNT_DAMAGE;
    }
    rfc_ctx->internal.flags                 = flags;

    /* Counter increments */
    rfc_ctx->full_inc                       = RFC_FULL_CYCLE_INCREMENT;
    rfc_ctx->half_inc                       = RFC_HALF_CYCLE_INCREMENT;
    rfc_ctx->curr_inc                       = RFC_FULL_CYCLE_INCREMENT;

    if( class_count )
    {
        if( class_count > RFC_CLASS_COUNT_MAX || class_width <= 0.0 )
        {
            return error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
    }
    else
    {
        class_width  = 1.0;
        class_offset = 0.0;
    }

    /* Rainflow class parameters */
    rfc_ctx->class_count                    = class_count;
    rfc_ctx->class_width                    = class_width;
    rfc_ctx->class_offset                   = class_offset;
    rfc_ctx->hysteresis                     = hysteresis;

    /* Values for a "pseudo Woehler curve" */
    rfc_ctx->state = RFC_STATE_INIT;   /* Bypass sanity check for state in wl_init() */
    RFC_wl_init_elementary( rfc_ctx, /*sx*/ RFC_WL_SD_DEFAULT, /*nx*/ RFC_WL_ND_DEFAULT, /*k*/ RFC_WL_K_DEFAULT );
    rfc_ctx->state = RFC_STATE_INIT0;  /* Reset state */

    /* Memory allocator */
    if( !rfc_ctx->mem_alloc )
    {
        rfc_ctx->mem_alloc = mem_alloc;
    }
    
    /* Residue */
    rfc_ctx->internal.residue_cap           = NUMEL( rfc_ctx->internal.residue );
    rfc_ctx->residue_cnt                    = 0;
    rfc_ctx->residue_cap                    = 2 * rfc_ctx->class_count + 1; /* 4pt-method fills max 2*n-2 (+1 candidate), +2 extra points ("enforce margin", "interim point") = 2*n+1 */

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


    rfc_ctx->state = RFC_STATE_INIT;

    return true;
}

/**
 * @brief      Return state
 *
 * @param      ctx   The rfc context
 *
 * @return     state
 */
rfc_state_e RFC_state_get( const void *ctx )
{
    /* RFC_CTX_CHECK_AND_ASSIGN */
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;                                           
                                                                                    
    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )                         
    {                                                                               
        return error_raise( rfc_ctx, RFC_ERROR_INVARG );                            
    }                                                                               

    return rfc_ctx->state;
}


/**
 * @brief      Return error
 *
 * @param      ctx   The rfc context
 *
 * @return     error
 */
rfc_error_e RFC_error_get( const void *ctx )
{
    RFC_CTX_CHECK_AND_ASSIGN

    return rfc_ctx->error;
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

/* Woehler curve */
    rfc_ctx->wl_sx        =  sx;
    rfc_ctx->wl_nx        =  nx;
    rfc_ctx->wl_k         = -fabs(k);

    return true;
}


/**
 * @brief      Returns the residuum
 *
 * @param      ctx              The rainflow context
 * @param[out] residue          The residue (last point is interim, if its tp_pos is zero)
 * @param[out] count            The residue count
 *
 * @return     true on success
 */
bool RFC_res_get( const void *ctx, const rfc_value_tuple_s **residue, unsigned *count )
{
    RFC_CTX_CHECK_AND_ASSIGN

    if( rfc_ctx->state < RFC_STATE_INIT )
    {
        return false;
    }

    if( residue )
    {
        *residue = rfc_ctx->residue;
    }

    if( count )
    {
        *count = (unsigned)rfc_ctx->residue_cnt + (rfc_ctx->state == RFC_STATE_BUSY_INTERIM);
    }

    return true;
}


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

    rfc_ctx->residue                    = NULL;
    rfc_ctx->residue_cap                = 0;
    rfc_ctx->residue_cnt                = 0;

    rfc_ctx->rfm                        = NULL;
    
    rfc_ctx->internal.slope             = 0;
    rfc_ctx->internal.extrema[0]        = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]        = nil;  /* local maximum */
    rfc_ctx->internal.pos               = 0;
    rfc_ctx->internal.pos_offset        = 0;

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

    if( !data ) return !data_count;

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

        if( rfc_ctx->class_count && ( tp.cls >= rfc_ctx->class_count || tp.value < rfc_ctx->class_offset ) )
        {
            return error_raise( rfc_ctx, RFC_ERROR_DATA_OUT_OF_RANGE );
        }
        
        if( !feed_once( rfc_ctx, &tp, rfc_ctx->internal.flags ) )
        {
            return false;
        }
    }

    return true;
}


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

    {
        int flags = rfc_ctx->internal.flags;

        switch( residual_method )
        {
            case RFC_RES_NONE:
                /* FALLTHROUGH */
            case RFC_RES_IGNORE:
                ok = finalize_res_ignore( rfc_ctx, flags );
                break;
            case RFC_RES_NO_FINALIZE:
                ok = finalize_res_no_finalize( rfc_ctx, flags );
                break;
            default:
                assert( false );
                ok = error_raise( rfc_ctx, RFC_ERROR_INVARG );
        }
        assert( rfc_ctx->state == RFC_STATE_FINALIZE );
    }

    if( !rfc_ctx->class_count )
    {
        rfc_ctx->residue_cnt = 0;
    }

    rfc_ctx->damage_residue = rfc_ctx->damage - damage;
    rfc_ctx->state          = ok ? RFC_STATE_FINISHED : RFC_STATE_ERROR;

    return ok;
}










/*** Implementation static functions ***/


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

    /* Check for next turning point and update residue. tp_residue is NULL, if there is no turning point */
    /* Otherwise tp_residue refers the forelast element in member rfc_ctx->residue */
    tp_residue = feed_filter_pt( rfc_ctx, pt );

    /* Countings */

    /* Add turning point and check for closed cycles */
    if( tp_residue )
    {

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

        if( tp_interim )
        {
            int flags = rfc_ctx->internal.flags;

            /* Check once more if a new cycle is closed now */
            cycle_find( rfc_ctx, flags );
        }

        rfc_ctx->state = RFC_STATE_FINALIZE;
    }

    return true;
}


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


/**
 * @brief      Finalize pending counts, ignore residue, don't finalize
 *
 * @param      rfc_ctx  The rainflow context
 * @param      flags    The flags
 *
 * @return     true on success
 */
static
bool finalize_res_no_finalize( rfc_ctx_s *rfc_ctx, rfc_flags_e flags )
{
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state < RFC_STATE_FINISHED );

    rfc_ctx->state = RFC_STATE_FINALIZE;

    return true;
}


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

            /* If D is integer format count:
             * 
             * D_integer = (double)(int)range ^ ABS(k)  // where range is 0..class_count-1
             * 
             * D = D_integer_sum * exp( log( ND * 2*SD/class_width ) * -ABS(k) )
             * 
             */

            /* Constants for the Woehler curve */
            const double SX_log = log(rfc_ctx->wl_sx);
            const double NX_log = log(rfc_ctx->wl_nx);
            const double k      = rfc_ctx->wl_k;

            /* Miner original */
            D = exp( fabs(k)  * ( log(Sa) - SX_log ) - NX_log );
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


/**
 * @brief      Calculate pseudo damage for one closed (full) cycle.
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
    double Sa = -1.0;  /* Negative amplitude states undefined */
    double D  =  0.0;

    assert( damage );
    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT );


    if( class_from != class_to )
    {
        Sa = fabs( (int)class_from - (int)class_to ) / 2.0 * rfc_ctx->class_width;

        if( !damage_calc_amplitude( rfc_ctx, Sa, &D ) )
        {
            return false;
        }
    }

    if( Sa_ret )
    {
        *Sa_ret = Sa;
    }

    *damage = D;

    return true;
}


/**
 * @brief      Test data sample for a new turning point and add to the residue
 *             in that case. Update extrema.
 *             - 1. Hysteresis Filtering
 *             - 2. Peak-Valley Filtering
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  pt       The data tuple, must not be NULL
 *
 * @return     Returns pointer to new turning point in residue or NULL
 */
static
rfc_value_tuple_s * feed_filter_pt( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *pt )
{
    int                 slope;
    rfc_value_t         delta;
    rfc_value_tuple_s  *new_tp      = NULL;
    bool                do_append   = false;

    assert( rfc_ctx );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state <= RFC_STATE_BUSY_INTERIM );

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
            }
            else if( pt->value > rfc_ctx->internal.extrema[1].value )
            {
                /* Maximum */
                is_falling_slope = 0;
                rfc_ctx->internal.extrema[1] = *pt;
            }

            /* Local hysteresis filtering */
            delta = value_delta( rfc_ctx, &rfc_ctx->internal.extrema[0], &rfc_ctx->internal.extrema[1], NULL /* sign_ptr */ );

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
                do_append = true;
            }
        }
    }
    else  /* if( rfc_ctx->state < RFC_STATE_BUSY_INTERIM ) */
    {
        assert( rfc_ctx->state == RFC_STATE_BUSY_INTERIM );

        /* Consecutive search for turning points */

        /* Hysteresis Filtering, check against interim turning point */
        delta = value_delta( rfc_ctx, &rfc_ctx->residue[rfc_ctx->residue_cnt], pt, &slope /* sign_ptr */ );

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
                do_append = true;
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
        assert( rfc_ctx->residue_cnt + 1 < rfc_ctx->residue_cap );
        rfc_ctx->residue[++rfc_ctx->residue_cnt] = *pt;

        /* Return new turning point */
        new_tp = &rfc_ctx->residue[rfc_ctx->residue_cnt - 1];
    }

    return new_tp;
}


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

    if( !rfc_ctx->class_count || ( from->value >= rfc_ctx->class_offset && to->value >= rfc_ctx->class_offset ) )
    {
        /* If class_count is zero, no counting is done. Otherwise values must be greater than class_offset */
    }
    else
    {
        assert( false );
    }

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

    }
}


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
 * @param      rfc_ctx   The rainflow context
 * @param[in]  pt_from   The point from
 * @param[in]  pt_to     The point to
 * @param[out] sign_ptr  Pointer to catch sign (may be NULL)
 *
 * @return     Returns the absolute difference of given values
 */
static
rfc_value_t value_delta( rfc_ctx_s* rfc_ctx, const rfc_value_tuple_s* pt_from, const rfc_value_tuple_s* pt_to, int *sign_ptr )
{
    double delta;

    assert( rfc_ctx );
    assert( pt_from && pt_to );

#if RFC_USE_HYSTERESIS_FILTER
    delta = (double)pt_to->value - (double)pt_from->value;
#else /*RFC_USE_HYSTERESIS_FILTER*/
    delta = rfc_ctx->class_width * ( (int)pt_to->cls - (int)pt_from->cls );
#endif /*RFC_USE_HYSTERESIS_FILTER*/

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
void * mem_alloc( void *ptr, size_t num, size_t size, int aim )
{
    if( !num || !size )
    {
        if( ptr )
        {
            FREE( ptr );
        }
        return NULL;
    }
    else
    {
        return ptr ? REALLOC( ptr, num * size ) : CALLOC( num, size );
    }
}
