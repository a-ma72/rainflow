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
 * This implementation uses the 4-point algorithm mentioned in [2] and [3].
 * To take the residue into account, you may implement a custom method or use some
 * predefined functions.
 * 
 * References:
 * [1] ASTM Standard E 1049, 1985 (2011). 
 *     "Standard Practices for Cycle Counting in Fatigue Analysis."
 *     West Conshohocken, PA: ASTM International, 2011.
 * [2] FVA-Richtlinie, 2010.
 *     "Zaehlverfahren zur Bildung von Kollektiven und Matrizen aus Zeitfunktionen"
 *     [https://fva-net.de/fileadmin/content/Richtlinien/FVA-Richtlinie_Zaehlverfahren_2010.pdf]
 * [3] Siemens Product Lifecycle Management Software Inc., 2018. 
 *     [https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Rainflow-Counting/ta-p/383093]
 * [4] G.Marsh on: "Review and application of Rainflow residue processing techniques for accurate fatigue damage estimation"
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

#include "rainflow.h"

#include <assert.h>  /* assert() */
#include <math.h>    /* exp(), log(), fabs() */


#ifdef MATLAB_MEX_FILE
#define RFC_MEX_USAGE \
"\nUsage:\n"\
"[pd,re,rm] = rfc( data, class_count, class_width, class_offset, hysteresis )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"
#pragma message(RFC_MEX_USAGE)
#include <string.h>
#include <mex.h>
#endif

/* Core functions */
static bool                 RFC_feed_once                       ( rfc_ctx_s *, rfc_value_tuple_s* tp );
static bool                 RFC_feed_finalize                   ( rfc_ctx_s *rfc_ctx );
static rfc_value_tuple_s *  RFC_tp_next                         ( rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 RFC_cycle_find_4ptm                 ( rfc_ctx_s * );
static void                 RFC_cycle_process                   ( rfc_ctx_s *, const rfc_value_tuple_s *from, const rfc_value_tuple_s *to, int flags );
/* Residual methods */
static bool                 RFC_finalize_res_ignore             ( rfc_ctx_s * );
/* Memory allocator */
static void *               RFC_mem_alloc                       ( void *ptr, size_t num, size_t size );
/* Other */
static double               RFC_damage_calc                     ( rfc_ctx_s *, unsigned class_from, unsigned class_to );
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
 *
 * @return     false on error
 */
bool RFC_init( void *ctx, 
               unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
               RFC_value_type hysteresis )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil   = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;
        return false;
    }

    if( rfc_ctx->state != RFC_STATE_INIT0 ) return false;

    /* Flags */
    rfc_ctx->flags                   = RFC_FLAGS_COUNT_MATRIX;
    
    /* Counter increments */
    rfc_ctx->full_inc                = RFC_FULL_CYCLE_INCREMENT;
    rfc_ctx->half_inc                = RFC_HALF_CYCLE_INCREMENT;
    rfc_ctx->curr_inc                = RFC_FULL_CYCLE_INCREMENT;

    if( !class_count || class_count > 512 ||
         class_width <= 0.0 )
    {
        rfc_ctx->error = RFC_ERROR_INVARG;
        return false;
    }

    /* Rainflow class parameters */
    rfc_ctx->class_count             = class_count;
    rfc_ctx->class_width             = class_width;
    rfc_ctx->class_offset            = class_offset;
    rfc_ctx->hysteresis              = ( hysteresis < 0.0 ) ? class_width : hysteresis;

    /* Woehler curve (fictive) */
    rfc_ctx->wl_sd                   =  1e3;            /* Fictive amplitude */
    rfc_ctx->wl_nd                   =  1e7;            /* Fictive count */
    rfc_ctx->wl_k                    = -5.0;            /* Fictive gradient */

    /* Memory allocator */
    if( !rfc_ctx->mem_alloc )
    {
        rfc_ctx->mem_alloc = RFC_mem_alloc;
    }

    /* Residue */
    rfc_ctx->residue_cnt             = 0;
    rfc_ctx->residue_cap             = 2 * rfc_ctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */
    rfc_ctx->residue                 = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, rfc_ctx->residue_cap, sizeof(rfc_value_tuple_s) );

    /* Non-sparse storages (optional, may be NULL) */
    rfc_ctx->matrix                  = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count * class_count, sizeof(RFC_counts_type) );

    /* Damage */
    rfc_ctx->pseudo_damage           = 0.0;

    rfc_ctx->internal.slope          = 0;
    rfc_ctx->internal.extrema[0]     = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]     = nil;  /* local maximum */

    if( !rfc_ctx->residue || !rfc_ctx->matrix )
    {
        RFC_deinit( rfc_ctx );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return false;
    }
    
    rfc_ctx->state = RFC_STATE_INIT;
    return true;
}


/**
 * @brief      De-initialization (freeing memory).
 *
 * @param      ctx  The rainflow context
 */
void RFC_deinit( void *ctx )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx || rfc_ctx->version != sizeof(rfc_ctx_s) )
    {
        assert( false );
        rfc_ctx->error = RFC_ERROR_INVARG;

        return;
    }

    if( rfc_ctx->residue )           rfc_ctx->mem_alloc( rfc_ctx->residue, 0, 0 );
    if( rfc_ctx->matrix )            rfc_ctx->mem_alloc( rfc_ctx->matrix,  0, 0 );

    rfc_ctx->residue                 = NULL;
    rfc_ctx->residue_cap             = 0;
    rfc_ctx->residue_cnt             = 0;

    rfc_ctx->matrix                  = NULL;

    rfc_ctx->internal.slope          = 0;
    rfc_ctx->internal.extrema[0]     = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]     = nil;  /* local maximum */
    rfc_ctx->internal.pos            = 0;


    rfc_ctx->state = RFC_STATE_INIT0;
}


/**
 * @brief      "Feed" counting algorithm with data samples (consecutive calls allowed).
 *
 * @param      ctx          The rainflow context
 * @param[in]  data         The data
 * @param[in]  data_count   The data count
 * 
 * @return     false on error
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

        if( !RFC_feed_once( rfc_ctx, &tp ) ) return false;
    }

    return true;
}


/**


/**
 * @brief       Finalize pending counts and turning point storage.
 *
 * @param       ctx              The rainflow context
 * @param       residual_method  The residual method (RFC_RES_...)
 * 
 * @return      false on error
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
    
    assert( rfc_ctx->state < RFC_STATE_FINALIZE );

    switch( residual_method )
    {
        case RFC_RES_NONE:
            /* fallthrough */
        case RFC_RES_IGNORE:
            ok = RFC_finalize_res_ignore( rfc_ctx );
            break;
        default:
            assert( false );
            rfc_ctx->error = RFC_ERROR_INVARG;
            ok = false;
    }
    
    assert( rfc_ctx->state == RFC_STATE_FINALIZE );

    rfc_ctx->state = ok ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
    return ok;
}







/*** Implementation static functions ***/

/**
 * @brief      Processing one data point.
 *             Find turning points and check for closed cycles.
 *
 * @param      rfc_ctx          The rainflow context
 * @param[in]  pt               The data tuple
 * 
 * @return     false on error
 */
static
bool RFC_feed_once( rfc_ctx_s *rfc_ctx, rfc_value_tuple_s* pt )
{
    rfc_value_tuple_s *tp_residue;

    assert( rfc_ctx && pt );

    /* Check for next turning point and update residue */
    /* (tp_residue refers member residue in rfc_ctx) */
    tp_residue = RFC_tp_next( rfc_ctx, pt );

    /* Rainflow counting */

    /* Check for closed cycles */
    if( tp_residue )
    {
        /* New turning point, check for closed cycles and count */
        RFC_cycle_find_4ptm( rfc_ctx );
    }

    return true;
}


/**
 * @brief      Handling interim turning point and margin.
 *             There are still unhandled turning point left.
 *             "Finalizing" takes this into account for the rainflow algorithm.
 *
 * @param      rfc_ctx  The rainflow context
 *
 * @return     false on error
 */
static
bool RFC_feed_finalize( rfc_ctx_s *rfc_ctx )
{
    rfc_value_tuple_s *tp_interim = NULL;

    assert( rfc_ctx );
    
    if( rfc_ctx->state >= RFC_STATE_FINALIZE )
    {
        return false;
    }

    /* Adjust residue: Incorporate interim turning point */
    if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
    {
        assert( rfc_ctx->residue && rfc_ctx->residue_cnt );

        tp_interim = &rfc_ctx->residue[rfc_ctx->residue_cnt];
        rfc_ctx->residue_cnt++;
    }

    if( tp_interim )
    {
        /* Check once more if a new cycle is closed now */
        RFC_cycle_find_4ptm( rfc_ctx );
    }

    rfc_ctx->state = RFC_STATE_FINALIZE;

    return true;
}


/**
 * @brief       Finalize pending counts, ignoring residue.
 *
 * @param       rfc_ctx  The rainflow context
 * 
 * @return      false on error
 */
static
bool RFC_finalize_res_ignore( rfc_ctx_s *rfc_ctx )
{
    /* Just include interim turning point */
    return RFC_feed_finalize( rfc_ctx );
}



/**
 * @brief      Calculate fictive damage for one closed (full) cycle.
 *
 * @param      rfc_ctx       The rainflow context
 * @param[in]  class_from    The starting class
 * @param[in]  class_to      The ending class
 *
 * @return     Pseudo damage value for the closed cycle
 */
static
double RFC_damage_calc( rfc_ctx_s *rfc_ctx, unsigned class_from, unsigned class_to )
{
    assert( rfc_ctx );


    /* Constants for the woehler curve */
    const double SD_log  = log(rfc_ctx->wl_sd);
    const double ND_log  = log(rfc_ctx->wl_nd);
    const double k       = rfc_ctx->wl_k;
    /* Pseudo damage */
    double D_i = 0.0;

    if( class_from != class_to )
    {
        /* D_i =           h_i /    ND   *    ( Sa_i /    SD)  ^ ABS(k)   */
        /* D_i = exp(  log(h_i /    ND)  + log( Sa_i /    SD)  * ABS(k) ) */
        /* D_i = exp( (log(h_i)-log(ND)) + (log(Sa_i)-log(SD)) * ABS(k) ) */
        /* D_i = exp(      0   -log(ND)  + (log(Sa_i)-log(SD)) * ABS(k) ) */

        double range  = (double)rfc_ctx->class_width * abs( (int)class_to - (int)class_from );
        double Sa_i   = range / 2.0;  /* amplitude */

        D_i = exp( fabs(k)  * ( log(Sa_i) - SD_log ) - ND_log );
    }

    return D_i;
}


/**
 * @brief       Returns the unsigned difference of two values, sign optionally returned as -1 or 1.
 *
 * @param[in]   from      Left hand value
 * @param[in]   to        Right hand value
 * @param[out]  sign_ptr  Pointer to catch sign (may be NULL)
 *
 * @return      Returns the absolute difference of given values
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
 * @brief      Test data sample for a new turning point and add to the residue in that case.
 *             1. Hysteresis Filtering
 *             2. Peak-Valley Filtering
 *
 * @param      rfc_ctx          The rainflow context
 * @param[in]  pt               The data point, must not be NULL
 *
 * @return     Returns pointer to new turning point in residue or NULL
 */
static
rfc_value_tuple_s * RFC_tp_next( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *pt )
{
    int                 slope;
    RFC_value_type      delta;
    rfc_value_tuple_s  *new_tp    = NULL;
    int                 do_append = 0;

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
            int is_falling_slope;

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
            delta = value_delta( rfc_ctx->internal.extrema[0].value, rfc_ctx->internal.extrema[1].value, NULL /* sign_ptr */ );

            if( delta > rfc_ctx->hysteresis )
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


/**
 * @brief      Rainflow counting core (4-point-method).
 *
 * @param      rfc_ctx  The rainflow context
 */
static
void RFC_cycle_find_4ptm( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx );

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

            RFC_cycle_process( rfc_ctx, from, to, rfc_ctx->flags );

            /* Remove two inner turning points (idx+1 and idx+2) */
            rfc_ctx->residue[idx+1] = rfc_ctx->residue[idx+3];  /* Move last turning point */
            rfc_ctx->residue[idx+2] = rfc_ctx->residue[idx+4];  /* Move interim turning point */
            rfc_ctx->residue_cnt -= 2;
        }
        else break;
    }
}


/**
 * @brief      Processes counts on a closing cycle
 *
 * @param      rfc_ctx  The rainflow context
 * @param[in]  from     The starting data point
 * @param[in]  to       The ending data point
 * @param[in]  flags    Control flags
 */
static
void RFC_cycle_process( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *from, const rfc_value_tuple_s *to, int flags )
{
    unsigned class_from, class_to;

    assert( rfc_ctx );
    assert( from->value > rfc_ctx->class_offset && to->value > rfc_ctx->class_offset );

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
        double damage = RFC_damage_calc( rfc_ctx, class_from, class_to );

        /* Adding damage due to current cycle weight */
        rfc_ctx->pseudo_damage += damage * rfc_ctx->curr_inc / rfc_ctx->full_inc;

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
    }
}


/**
 * @brief       (Re-)Allocate or free memory
 * @param       ptr     Previous data pointer, or NULL, if unset
 * @param[in]   num     The number of elements
 * @param[in]   size    The size of one element in bytes
 *
 * @returns     New memory pointer or NULL if either num or size is 0
 * 
 */
static
void * RFC_mem_alloc( void *ptr, size_t num, size_t size )
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




#ifdef MATLAB_MEX_FILE
/**
 * MATLAB wrapper for the rainflow algorithm
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    if( !nrhs )
    {
        mexPrintf( "%s", RFC_MEX_USAGE );
        return;
    }
    
    if( nrhs != 5 )
    {
        mexErrMsgTxt( "Function needs exact 5 arguments!" );
    }
    else
    {
        rfc_ctx_s rfc_ctx = { sizeof(rfc_ctx_s) };
    
        const mxArray *mxData        = prhs[0];
        const mxArray *mxClassCount  = prhs[1];
        const mxArray *mxClassWidth  = prhs[2];
        const mxArray *mxClassOffset = prhs[3];
        const mxArray *mxHysteresis  = prhs[4];

        RFC_value_type *buffer       = NULL;
        double         *data         = mxGetPr( mxData );
        size_t          data_len     = mxGetNumberOfElements( mxData );
        unsigned        class_count  = (unsigned)( mxGetScalar( mxClassCount ) + 0.5 );
        double          class_width  = mxGetScalar( mxClassWidth );
        double          class_offset = mxGetScalar( mxClassOffset );
        double          hysteresis   = mxGetScalar( mxHysteresis );
        size_t          i;
        bool            ok;

        ok = RFC_init( &rfc_ctx, 
                       class_count, (RFC_value_type)class_width, (RFC_value_type)class_offset, 
                       (RFC_value_type)hysteresis );

        if( !ok )
        {
            mexErrMsgTxt( "Error during initialization!" );
        }

        /* Casting values from double type to RFC_value_type */ 
        if( sizeof( RFC_value_type ) != sizeof(double) && data_len )  /* maybe unsafe! */
        {
            buffer = (RFC_value_type *)RFC_mem_alloc( NULL, data_len, sizeof(RFC_value_type) );

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
        RFC_feed( &rfc_ctx, buffer, data_len  );
        RFC_finalize( &rfc_ctx, RFC_RES_IGNORE );

        /* Free temporary buffer (cast) */
        if( (void*)buffer != (void*)data )
        {
            RFC_mem_alloc( buffer, 0, 0 );
            buffer = NULL;
        }

        /* Return results */
        if( plhs )
        {
            /* Pseudo damage */
            plhs[0] = mxCreateDoubleScalar( rfc_ctx.pseudo_damage );

            /* Residue */
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

            /* Rainflow matrix (column major order) */
            if( nlhs > 2 && rfc_ctx.matrix )
            {
                mxArray* matrix = mxCreateDoubleMatrix( class_count, class_count, mxREAL );
                if( matrix )
                {
                    mxArray* transposed = NULL;

                    if( sizeof( RFC_counts_type ) == sizeof(double) )  /* maybe unsafe! */
                    {
                        memcpy( mxGetPr(matrix), rfc_ctx.matrix, sizeof(double) * class_count * class_count );
                        mexCallMATLAB( 1, &transposed, 1, &matrix, "transpose" );
                        mxDestroyArray( matrix );
                    }
                    else
                    {
                        double *ptr = mxGetPr(matrix);
                        size_t from, to;
                        for( to = 0; to < class_count; to++ )
                        {
                            for( from = 0; from < class_count; from++ )
                            {
                                *ptr++ = (double)rfc_ctx.matrix[ from * class_count + to ];
                            }
                        }
                        transposed = matrix;
                    }

                    if( transposed )
                    {
                        plhs[2] = transposed;
                    }
                }
            }
        }

        /* Deinitialize rainflow context */
        RFC_deinit( &rfc_ctx );
    }
}
#endif
