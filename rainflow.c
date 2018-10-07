/** Rainflow Counting Algorithm (4-point-method), C99 compliant
 * 
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
"[pd,re,rm,rp,lc,tp] = rfc( data, class_count, class_width, class_offset, hysteresis )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"\
"    rp = Range pair counts\n"\
"    lc = Level crossings\n"\
"    tp = Turning points\n"
#pragma message(RFC_MEX_USAGE)
#include <string.h>
#include <mex.h>
#endif

/* Core */
static void             RFC_feed_handle_tp                  ( rfctx_s *rfctx, value_tuple_s* tp, bool is_last );
static void             RFC_feed_finalize                   ( rfctx_s* rfctx );
static value_tuple_s *  RFC_tp_next_default                 ( rfctx_s *, value_tuple_s *pt );
static void             RFC_cycle_find_4ptm                 ( rfctx_s * );
static void             RFC_cycle_process                   ( rfctx_s *, value_tuple_s *from, value_tuple_s *to, int flags );
/* Residual methods */
static bool             RFC_finalize_default                ( rfctx_s * );
static bool             RFC_finalize_res_default            ( rfctx_s * );
static bool             RFC_finalize_res_ignore             ( rfctx_s * );
static bool             RFC_finalize_res_halfcycles         ( rfctx_s * );
static bool             RFC_finalize_res_fullcycles         ( rfctx_s * );
static bool             RFC_finalize_res_clormann_seeger    ( rfctx_s * );
static bool             RFC_finalize_res_din                ( rfctx_s * );
static bool             RFC_finalize_res_repeated           ( rfctx_s * );
/* Other */
static void             RFC_tp_add_default                  ( rfctx_s *, value_tuple_s *pt, bool do_lock );
static double           RFC_damage_calc_default             ( rfctx_s *, unsigned class_from, unsigned class_to );
static RFC_value_type   value_delta                         ( RFC_value_type from, RFC_value_type to, int *sign_ptr );


#define QUANTIZE( r, v )   ( (unsigned)( ((v) - (r)->class_offset) / (r)->class_width ) )
#define CLASS_MEAN( r, c ) ( (double)( (r)->class_width * (0.5 + (c)) + (r)->class_offset ) )


/**
 * @brief      Initialization (rainflow context).
 *
 * @param      ctx           The rainflow context
 * @param[in]  class_count   The class count
 * @param[in]  class_width   The class width
 * @param[in]  class_offset  The class offset
 * @param[in]  hysteresis    The hysteresis
 * @param[in]  tp            Pointer to turning points buffer
 * @param[in]  tp_cap        Number of turning points in buffer
 *
 * @return     true on success
 */
bool RFC_init( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                          RFC_value_type hysteresis, 
                          int residual_method,
                          value_tuple_s *tp, size_t tp_cap )
{
    rfctx_s       *rfctx = (rfctx_s*)ctx;
    value_tuple_s  nil   = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfctx || rfctx->state != RFC_STATE_INIT0 ) return false;

    assert( rfctx->version == sizeof(rfctx_s) );

    /* Flags */
    rfctx->flags = RFC_FLAGS_COUNT_ALL;
    
    /* Residual method */
    rfctx->residual_method      = residual_method;

    /* Counter increments */
    rfctx->full_inc             = RFC_FULL_CYCLE_INCREMENT;
    rfctx->half_inc             = RFC_HALF_CYCLE_INCREMENT;
    rfctx->curr_inc             = RFC_FULL_CYCLE_INCREMENT;

    if( !class_count || class_count > 512 ) return false;
    if( class_width <= 0.0 )                return false;

    /* Rainflow class parameters */
    rfctx->class_count          = class_count;
    rfctx->class_width          = class_width;
    rfctx->class_offset         = class_offset;
    rfctx->hysteresis           = ( hysteresis < 0.0 ) ? class_width : hysteresis;

    /* Woehler curve (fictive) */
    rfctx->wl_sd                =  1e3;          /* Fictive value */
    rfctx->wl_nd                =  1e7;          /* Fictive value */
    rfctx->wl_k                 = -5.0;          /* Fictive value */
    rfctx->wl_k2                =  rfctx->wl_k;  /* Miner elementar, if k == k2 */
    rfctx->wl_omission          =  0.0;          /* No omission per default */

    /* Memory allocation functions */
    rfctx->mem_alloc            = calloc;
    rfctx->mem_free             = free;

#if RFC_USE_DELEGATES
    /* Delegates (optional, may be NULL) */
    rfctx->tp_next_fcn          = RFC_tp_next_default;
    rfctx->tp_add_fcn           = RFC_tp_add_default;
    rfctx->cycle_find_fcn       = RFC_cycle_find_4ptm;
    rfctx->finalize_fcn         = RFC_finalize_res_ignore;
    rfctx->damage_calc_fcn      = RFC_damage_calc_default;
#endif

    /* Residue */
    rfctx->residue_cnt          = 0;
    rfctx->residue_cap          = 2 * rfctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */
    rfctx->residue              = (value_tuple_s*)rfctx->mem_alloc( rfctx->residue_cap, sizeof(value_tuple_s) );

    /* Non-sparse storages (optional, may be NULL) */
    rfctx->matrix               = (RFC_counts_type*)rfctx->mem_alloc( class_count * class_count, sizeof(RFC_value_type) );
    rfctx->rp                   = (RFC_counts_type*)rfctx->mem_alloc( class_count,               sizeof(RFC_value_type) );
    rfctx->lc                   = (RFC_counts_type*)rfctx->mem_alloc( class_count,               sizeof(RFC_value_type) );

    /* Damage */
    rfctx->pseudo_damage        = 0.0;

    rfctx->internal.slope       = 0;
    rfctx->internal.extrema[0]  = nil;
    rfctx->internal.extrema[1]  = nil;
    rfctx->internal.margin[0]   = nil;
    rfctx->internal.margin[1]   = nil;

    if( !rfctx->residue || !rfctx->matrix || !rfctx->rp || !rfctx->lc )
    {
        RFC_deinit( rfctx );
        return false;
    }
    
    /* Turning points storage (optional, may be NULL) */
    rfctx->tp                   = tp;
    rfctx->tp_cap               = tp ? tp_cap : 0;
    rfctx->tp_cnt               = 0;
    rfctx->tp_locked            = 0;

    rfctx->state = RFC_STATE_INIT;
    return true;
}


/**
 * @brief      De-initialization (freeing memory).
 *
 * @param      ctx  The rainflow context
 */
void RFC_deinit( void *ctx )
{
    rfctx_s       *rfctx = (rfctx_s*)ctx;
    value_tuple_s  nil   = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfctx ) return;

    assert( rfctx->version == sizeof( rfctx_s ) );

    if( rfctx->residue )       rfctx->mem_free( rfctx->residue );
    if( rfctx->matrix )        rfctx->mem_free( rfctx->matrix );
    if( rfctx->rp )            rfctx->mem_free( rfctx->rp );
    if( rfctx->lc )            rfctx->mem_free( rfctx->lc );
    if( rfctx->tp )            rfctx->mem_free( rfctx->tp );

    rfctx->residue             = NULL;
    rfctx->residue_cap         = 0;
    rfctx->residue_cnt         = 0;

    rfctx->matrix              = NULL;
    rfctx->rp                  = NULL;
    rfctx->lc                  = NULL;

    rfctx->internal.slope      = 0;
    rfctx->internal.extrema[0] = nil;
    rfctx->internal.extrema[1] = nil;
    rfctx->internal.margin[0]  = nil;
    rfctx->internal.margin[1]  = nil;
    rfctx->internal.pos        = 0;

    rfctx->tp                  = NULL;
    rfctx->tp_cap              = 0;
    rfctx->tp_cnt              = 0;
    rfctx->tp_locked           = 0;

    rfctx->state = RFC_STATE_INIT0;
}


/**
 * @brief      Feed counting algorithm with data samples.
 *
 * @param      ctx          The rainflow context
 * @param[in]  data         The data
 * @param[in]  data_count   The data count
 * @param[in]  do_finalize  Flag to finalize counting
 */
void RFC_feed( void *ctx, const RFC_value_type * data, size_t data_count, bool do_finalize )
{
    rfctx_s *rfctx = (rfctx_s*)ctx;

    if( !rfctx ) return;

    assert( !data_count || data );
    assert( rfctx->version == sizeof( rfctx_s ) );

    if( rfctx->state < RFC_STATE_INIT || rfctx->state >= RFC_STATE_FINISHED )
    {
        assert( false );
        return;
    }

    /* Process data */
    while( data_count-- )
    {
        value_tuple_s tp = { (RFC_value_type)*data++ };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

        /* Assign class and global position (base 1) */
        tp.class = QUANTIZE( rfctx, tp.value );
        tp.pos   = ++rfctx->internal.pos;

        RFC_feed_handle_tp( rfctx, &tp, do_finalize && !data_count /* is_last */ );
    }

    /* Last data point? */
    if( do_finalize )
    {
        RFC_feed_finalize( rfctx );
    }
}


/**
 * @brief      Feed counting algorithm with data tuples.
 *
 * @param      ctx          The rainflow context
 * @param[in]  data         The data tuples
 * @param[in]  data_count   The data count
 * @param[in]  do_finalize  Flag to finalize counting
 */
void RFC_feed_tuple( void *ctx, value_tuple_s *data, size_t data_count, bool do_finalize )
{
    rfctx_s *rfctx = (rfctx_s*)ctx;

    if( !rfctx ) return;

    assert( !data_count || data );
    assert( rfctx->version == sizeof( rfctx_s ) );

    if( rfctx->state < RFC_STATE_INIT || rfctx->state >= RFC_STATE_FINISHED )
    {
        assert( false );
        return;
    }

    /* Process data */
    while( data_count-- )
    {
        RFC_feed_handle_tp( rfctx, data++, do_finalize && !data_count );
    }

    /* Last data point? */
    if( do_finalize )
    {
        RFC_feed_finalize( rfctx );
    }
}


/**
 * @brief      Processing data (one data point).
 *             Find turning points and check for closed cycles.
 *
 * @param      rfctx        The rainflow context
 * @param[in]  tp           The data tuple
 * @param[in]  is_last      True, if this is the last data tuple
 */
static
void RFC_feed_handle_tp( rfctx_s *rfctx, value_tuple_s* tp, bool is_last )
{
    value_tuple_s *tp_residue;

    if( is_last )
    {
        rfctx->margin[1] = *tp;
    }

#if RFC_USE_DELEGATES
    /* Test if a new turning point exists */
    if( rfctx->tp_next_fcn && ( tp_residue = rfctx->tp_next_fcn( rfctx, tp ) ) )
    {
        if( rfctx->tp_add_fcn )
        {
            /* Add new turning point */
            rfctx->tp_add_fcn( rfctx, tp_residue, false /* do_lock */ );
        }

        /* New turning point, check for closed cycles and count */
        if( rfctx->cycle_find_fcn )
        {
            rfctx->cycle_find_fcn( rfctx );
        }
    }
#else
    /* Test if a new turning point exists */
    if( tp_residue = RFC_tp_next_default( rfctx, tp ) )
    {
        /* Add new turning point */
        RFC_tp_add_default( rfctx, tp_residue, false /* do_lock */ );

        /* New turning point, check for closed cycles and count */
        RFC_cycle_find_4ptm( rfctx );
    }
#endif
}


/**
 * @brief      Finalizing.
 *             Take residue into account.
 *
 * @param      rfctx        The rainflow context
 */
static
void RFC_feed_finalize( rfctx_s* rfctx )
{
    rfctx->state = RFC_STATE_FINALIZE;

    /* Adjust residue: Incorporate interim turning point */
    rfctx->residue_cnt++;

    /* Finalizing (process last turning point */
#if RFC_USE_DELEGATES
    if( rfctx->finalize_fcn )
    {
        rfctx->state = rfctx->finalize_fcn( rfctx ) ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
    }
#else
    rfctx->state = RFC_finalize_res_default( rfctx ) ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
#endif
}


/**
 * @brief      Finalize pending counts, handling interim turning point.
 *             There is still one unhandled turning point left.
 *             "Finalizing" takes this into account.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_default( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
#if RFC_USE_DELEGATES
        /* Adjust turning points: Append very last turning point and lock */
        /* Add new turning point */
        if( rfctx->tp_add_fcn )
        {
            /* Add new turning point */
            rfctx->tp_add_fcn( rfctx, &rfctx->residue[rfctx->residue_cnt - 1], true /* lock */ );
        }

        if( rfctx->cycle_find_fcn )
        {
            /* Check if a new cycle is closed now */
            rfctx->cycle_find_fcn( rfctx );
        }
#else
        /* Add interim turning point */
        RFC_tp_add_default( rfctx, &rfctx->residue[rfctx->residue_cnt - 1], true /* lock */ );

        /* Check if a new cycle is closed now */
        RFC_cycle_find_4ptm( rfctx );
#endif
    }
    return true;
}


/**
 * @brief      Finalize pending counts, ignoring residue.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_default( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    switch( rfctx->residual_method )
    {
        case RFC_RES_NONE:
        case RFC_RES_IGNORE:
            return RFC_finalize_res_ignore( rfctx );
        case RFC_RES_HALFCYCLES:
            return RFC_finalize_res_halfcycles( rfctx );
        case RFC_RES_FULLCYCLES:
            return RFC_finalize_res_fullcycles( rfctx );
        case RFC_RES_CLORMANN_SEEGER:
            return RFC_finalize_res_clormann_seeger( rfctx );
        case RFC_RES_REPEATED:
            return RFC_finalize_res_repeated( rfctx );
        case RFC_RES_RP_DIN:
            return RFC_finalize_res_din( rfctx );
        default:
            assert( false );
            return false;
    }
}


/**
 * @brief      Finalize pending counts, ignoring residue.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_ignore( rfctx_s *rfctx )
{
    /* Just include interim turning point */
    RFC_finalize_default( rfctx );
}


/**
 * @brief      Finalize pending counts, half cycles method.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_halfcycles( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
    }
}


/**
 * @brief      Finalize pending counts, full cycles method.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_fullcycles( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
    }
}


/**
 * @brief      Finalize pending counts, Clormann/Seeger method.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_clormann_seeger( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
    }
}


/**
 * @brief      Finalize pending counts, DIN method.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_din( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
    }
}


/**
 * @brief      Finalize pending counts, repeated residue method.
 *
 * @param      rfctx  The rainflow context
 */
static
bool RFC_finalize_res_repeated( rfctx_s *rfctx )
{
    assert( rfctx && rfctx->state == RFC_STATE_FINALIZE );

    if( rfctx->residue && rfctx->residue_cnt )
    {
        /* Include interim turning point as new data series, 
           but don't modify residue history itself */
        size_t          cnt     = rfctx->residue_cnt;
        value_tuple_s  *residue = (value_tuple_s*)rfctx->mem_alloc( ++cnt, sizeof(value_tuple_s) );

        if( residue )
        {
            /* Make a copy of the residue */
            size_t n = cnt;
            const value_tuple_s *from = rfctx->residue;
                  value_tuple_s *to   = residue;

            while( n-- )
            {
                *to++ = *from++;
            }

            /* Make last turning point interim again */
            rfctx->residue_cnt--;
            
            /* Feed again with the copy */
            RFC_feed_tuple( rfctx, residue, cnt, false /*do_finalize*/ );

            /* Include interim turning point again */
            rfctx->residue_cnt++;

            rfctx->mem_free( residue );
        }
        else return false;

        /* Handle interim turning point */
        return RFC_finalize_default( rfctx );
    }
    return true;
}


/**
 * @brief      Calculate fictive damage for one closed (full) cycle.
 *
 * @param      rfctx       The rainflow context
 * @param[in]  class_from  The starting class
 * @param[in]  class_to    The ending class
 *
 * @return     Pseudo damage value for the closed cycle
 */
static
double RFC_damage_calc_default( rfctx_s *rfctx, unsigned class_from, unsigned class_to )
{
    assert( rfctx );

    /* Constants for woehler curve */
    const double SD_log  = log(rfctx->wl_sd);
    const double ND_log  = log(rfctx->wl_nd);
    const double k       = rfctx->wl_k;
    const double k2      = rfctx->wl_k2;
    /* Pseudo damage */
    double D_i = 0.0;

    if( class_from != class_to )
    {
        /* D_i =           h_i /    ND   *    ( Sa_i /    SD)  ^ ABS(k)   */
        /* D_i = exp(  log(h_i /    ND)  + log( Sa_i /    SD)  * ABS(k) ) */
        /* D_i = exp( (log(h_i)-log(ND)) + (log(Sa_i)-log(SD)) * ABS(k) ) */
        /* D_i = exp(      0   -log(ND)  + (log(Sa_i)-log(SD)) * ABS(k) ) */

        double range  = (double)rfctx->class_width * abs( (int)class_to - (int)class_from );
        double Sa_i   = range / 2.0;  /* amplitude */

        if( Sa_i > rfctx->wl_omission )
        {
            if( Sa_i > rfctx->wl_sd )
            {
                D_i = exp( fabs(k)  * ( log(Sa_i) - SD_log ) - ND_log );
            }
            else
            {
                D_i = exp( fabs(k2) * ( log(Sa_i) - SD_log ) - ND_log );
            }
        }
    }

    return D_i;
}


/*** Implementation static functions ***/

/**
 * @brief      Returns the unsigned difference of two values, sign optionally returned as -1 or 1.
 *
 * @param[in]  from      Left hand value
 * @param[in]  to        Right hand value
 * @param      sign_ptr  Pointer to catch sign (may be NULL)
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
 * @brief      Test data sample for a new turning point and add to the residue if so.
 *             1. Hysteresis Filtering
 *             2. Peak-Valley Filtering
 *
 * @param      rfctx  The rainflow context
 * @param      pt     The data point
 *
 * @return     Returns pointer to new turning point in residue or NULL
 */
static
value_tuple_s * RFC_tp_next_default( rfctx_s *rfctx, value_tuple_s *pt )
{
    int             slope;
    RFC_value_type  delta;

    assert( rfctx );

    slope = rfctx->internal.slope;

    if( !rfctx->residue_cnt )
    {
        /* Residue is empty, still searching the first turning point(s) */

        int i_first_point;  /* Index to internal.extrema array indicating 1st point */

        if( rfctx->state == RFC_STATE_INIT )
        {
            /* Very first point, initialize local min-max search and margin */
            rfctx->internal.extrema[0] = 
            rfctx->internal.extrema[1] =
            rfctx->internal.margin[0]  = *pt;
            rfctx->state               = RFC_STATE_BUSY;

            if( rfctx->flags & RFC_FLAGS_ENFORCE_MARGIN )
            {
                /* Enforce margin at the beginning */
                assert( rfctx->residue_cnt < rfctx->residue_cap );
                rfctx->residue[rfctx->residue_cnt++] = rfctx->internal.margin[0];
            }
        }
        else
        {
            /* Still searching for first turning point(s) */

            /* Update local extrema */
            if( pt->value < rfctx->internal.extrema[0].value )
            {
                /* Minimum */
                i_first_point = 1;
                rfctx->internal.extrema[0] = *pt;
            }
            else if( pt->value > rfctx->internal.extrema[1].value )
            {
                /* Maximum */
                i_first_point = 0;
                rfctx->internal.extrema[1] = *pt;
            }
        }

        /* Local hysteresis filtering */
        delta = value_delta( rfctx->internal.extrema[0].value, rfctx->internal.extrema[1].value, NULL /* sign_ptr */ );

        if( delta > rfctx->hysteresis )
        {
            /* Criteria met, first turning point found */
            assert( rfctx->residue_cnt < rfctx->residue_cap );
            rfctx->residue[rfctx->residue_cnt++] = rfctx->internal.extrema[i_first_point];

            /* Append interim turning point */
            assert( rfctx->residue_cnt < rfctx->residue_cap );
            rfctx->residue[rfctx->residue_cnt] = *pt;

            rfctx->internal.slope = i_first_point ? -1 : 1;

            return rfctx->residue;
        }
        
        return NULL;
    }

    /* Consecutive search for turning points */

    /* Hysteresis Filtering, check against interim turning point */
    delta = value_delta( rfctx->residue[rfctx->residue_cnt].value, pt->value, &slope /* sign_ptr */ );

    /* There are three scenarios possible here:
     *   1. Slope reversal, slope is less than or equal hysteresis
     *      Nothing to do.
     *   2. Previous slope is continued
     *      "delta" is ignored whilst hysteresis is exceeded already.
     *      Interim turning point has just to be adjusted.
     *   3. Slope reversal, slope is greater than hysteresis
     *      Interim turning point becomes real turning point.
     *      Current point becomes new interim turning point
     */

    /* Peak-Valley Filtering */
    if( slope != rfctx->internal.slope && delta <= rfctx->hysteresis )
    {
        /* Scenario (1), Turning point found, but still in hysteresis band */
        return NULL;
    }

    /* Justify interim turning point, or add a new one (2) */
    if( slope == rfctx->internal.slope )
    {
        /* Scenario (2), Continuous slope */

        /* Replace interim turning point with new extrema */
        rfctx->residue[rfctx->residue_cnt] = *pt;

        /* No new turning point */
        return NULL;
    }
    else
    {
        /* Scenario (3), Criteria met: slope != rfctx->internal.slope && delta > rfctx->hysteresis */

        /* New interim turning point */
        assert( rfctx->residue_cnt < rfctx->residue_cap );
        rfctx->residue[++rfctx->residue_cnt] = *pt;

        rfctx->internal.slope = slope;

        /* Have new turning point */
        return &rfctx->residue[rfctx->residue_cnt - 1];
    }
}


/**
 * @brief      Rainflow counting core (4-point-method).
 *
 * @param      rfctx  The rainflow context
 */
static
void RFC_cycle_find_4ptm( rfctx_s *rfctx )
{
    assert( rfctx );

    while( rfctx->residue_cnt >= 4 )
    {
        size_t idx = rfctx->residue_cnt - 4;

        RFC_value_type A = rfctx->residue[idx+0].value;
        RFC_value_type B = rfctx->residue[idx+1].value;
        RFC_value_type C = rfctx->residue[idx+2].value;
        RFC_value_type D = rfctx->residue[idx+3].value;

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
            value_tuple_s *from = &rfctx->residue[idx+1];
            value_tuple_s *to   = &rfctx->residue[idx+2];

            RFC_cycle_process( rfctx, from, to, rfctx->flags );

            /* Remove two inner turning points (idx+1 and idx+2)*/
            rfctx->residue[idx+1] = rfctx->residue[idx+3];  /* Move last turning point */
            rfctx->residue[idx+2] = rfctx->residue[idx+4];  /* Move interim turning point */
            rfctx->residue_cnt -= 2;
        }
        else break;
    }
}


/**
 * @brief      Processes counts on a closing cycle
 *
 * @param      rfctx  The rainflow context
 * @param      from   The starting data point
 * @param      to     The ending data point
 */
static
void RFC_cycle_process( rfctx_s *rfctx, value_tuple_s *from, value_tuple_s *to, int flags )
{
    unsigned class_from, class_to;

    assert( rfctx );
    assert( from->value > rfctx->class_offset && to->value > rfctx->class_offset );

    /* If flag RFC_FLAGS_ENFORCE_MARGIN is set, cycles less than hysteresis are possible */
    if( flags & RFC_FLAGS_ENFORCE_MARGIN )
    {
        if( value_delta( from->value, to->value, NULL /* sign_ptr */ ) <= rfctx->hysteresis )
        {
            return;
        }
    }

    /* Quantize "from" */
    class_from = QUANTIZE( rfctx, from->value );

    if( class_from >= rfctx->class_count ) class_from = rfctx->class_count - 1;

    /* Quantize "to" */
    class_to = QUANTIZE( rfctx, to->value );

    if( class_to >= rfctx->class_count ) class_to = rfctx->class_count - 1;
    
    /* class_from and class_to are base 0 now */

    /* Do several countings, according to "flags" */
    if( class_from != class_to )
    {
        /* Cumulate pseudo damage */
        double damage;
#if RFC_USE_DELEGATES
        if( rfctx->damage_calc_fcn )
        {
             damage = rfctx->damage_calc_fcn( rfctx, class_from, class_to );
        }
#else
        damage = RFC_damage_calc_default( rfctx, class_from, class_to );
#endif
        /* Adding damage due to current cycle weight */
        rfctx->pseudo_damage += damage * rfctx->curr_inc / rfctx->full_inc;

        /* Rainflow matrix */
        if( rfctx->matrix && ( flags & RFC_FLAGS_COUNT_MATRIX ) )
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
            size_t idx = rfctx->class_count * class_from + class_to;
            
            assert( rfctx->matrix[idx] <= RFC_COUNTS_LIMIT );
            rfctx->matrix[idx] += rfctx->curr_inc;
        }

        /* Range pair */
        if( rfctx->rp && ( flags & RFC_FLAGS_COUNT_RP ) )
        {
            /* 
             * Range pair histogram (vector storage)
             * Range value = idx * class_width  (=2x Amplitude)
             */
            int idx = abs( (int)class_from - (int)class_to );
            
            assert( rfctx->rp[idx] <= RFC_COUNTS_LIMIT );
            rfctx->rp[idx] += rfctx->curr_inc;
        }

        /* Level crossing, count rising and falling slopes */
        if( rfctx->lc && ( flags & RFC_FLAGS_COUNT_LC ) )
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
                    assert( rfctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfctx->lc[idx] += rfctx->curr_inc;
                }
            }
            else if( class_to < class_from && ( flags & RFC_FLAGS_COUNT_LC_DN ) )
            {
                /* Count falling slopes */
                unsigned idx;
                for( idx = class_to; idx < class_from; idx++ )
                {
                    assert( rfctx->lc[idx] <= RFC_COUNTS_LIMIT );
                    rfctx->lc[idx] += rfctx->curr_inc;
                }
            }
        }
    }
}


/**
 * Append one data sample to the turning points queue
 */
static
void RFC_tp_add_default( rfctx_s *rfctx, value_tuple_s *pt, bool do_lock )
{
    assert( rfctx );

    if( rfctx->tp && !rfctx->tp_locked )
    {
        if( pt )
        {
            assert( rfctx->tp_cnt < rfctx->tp_cap );
            rfctx->tp[ rfctx->tp_cnt++ ] = *pt;
        }

        if( do_lock )
        {
            rfctx->tp_locked = 1;
        }
    }
}


/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


void RFC_lc_from_matrix( rfctx_s *rfctx, RFC_counts_type* buffer, size_t buffer_len );
void RFC_rp_from_matrix( rfctx_s *rfctx, RFC_counts_type* buffer, size_t buffer_len );

/**
 * Calculate level crossing counts from rainflow matrix, write results to buffer.
 */
void RFC_lc_from_matrix( rfctx_s *rfctx, RFC_counts_type* buffer, size_t buffer_len )
{
    unsigned i, j, k;
    bool     up     = rfctx->flags & RFC_FLAGS_COUNT_LC_UP;
    bool     dn     = rfctx->flags & RFC_FLAGS_COUNT_LC_DN;
    size_t   maxcnt = buffer_len / sizeof(RFC_counts_type);

    assert( rfctx );

    if( !buffer || !maxcnt ) return;

    for( i = 0; i < rfctx->class_count; i++ ) 
    {
        RFC_counts_type counts = (RFC_counts_type)0;

        for( j = i; j < rfctx->class_count; j++ )  /* To */
        {
            for( k = 0; k < i; k++ )               /* From */
            {
                /* Count rising slopes */
                assert( counts < RFC_COUNTS_LIMIT - rfctx->matrix[ k * rfctx->class_count + j ] );
                if( up )
                {
                    counts += rfctx->matrix[ k * rfctx->class_count + j ];
                }

                /* Count falling slopes */
                assert( counts < RFC_COUNTS_LIMIT - rfctx->matrix[ j * rfctx->class_count + k ] );
                if( dn )
                {
                    counts += rfctx->matrix[ j * rfctx->class_count + k ];
                }
            }
        }

        buffer[i > maxcnt ? maxcnt : i] = counts;
    }
}


/**
 * Calculate range pair counts from rainflow matrix, write results to buffer.
 */
void RFC_rp_from_matrix( rfctx_s *rfctx, RFC_counts_type* buffer, size_t buffer_len )
{
    unsigned i, j;
    size_t   maxcnt = buffer_len / sizeof(RFC_counts_type);

    assert( rfctx && buffer && buffer_len >= rfctx->class_count * sizeof(RFC_counts_type) );

    buffer[0] = 0;

    for( i = 1; i < rfctx->class_count; i++ ) 
    {
        RFC_counts_type counts = (RFC_counts_type)0;

        for( j = i; j < rfctx->class_count; j++ )
        {
            /* Count rising and falling slopes */
            assert( counts < RFC_COUNTS_LIMIT - rfctx->matrix[ i * rfctx->class_count + j ] 
                                              - rfctx->matrix[ j * rfctx->class_count + i ] );
            counts += rfctx->matrix[ i * rfctx->class_count + j ];
            counts += rfctx->matrix[ j * rfctx->class_count + i ];
        }

        buffer[i > maxcnt ? maxcnt : i] = counts;
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
        rfctx_s rfctx = { sizeof(rfctx_s) };
    
        const mxArray *mxData        = prhs[0];
        const mxArray *mxClassCount  = prhs[1];
        const mxArray *mxClassWidth  = prhs[2];
        const mxArray *mxClassOffset = prhs[3];
        const mxArray *mxHysteresis  = prhs[4];

        double   *data         = mxGetPr( mxData );
        size_t    data_len     = mxGetNumberOfElements( mxData );
        unsigned  class_count  = (unsigned)( mxGetScalar( mxClassCount ) + 0.5 );
        double    class_width  = mxGetScalar( mxClassWidth );
        double    class_offset = mxGetScalar( mxClassOffset );
        double    hysteresis   = mxGetScalar( mxHysteresis );

        rfctx.tp = (value_tuple_s*)calloc( data_len, sizeof(value_tuple_s) );

        if( !RFC_init( &rfctx, class_count, class_width, class_offset, hysteresis, RFC_RES_NONE, rfctx.tp, data_len ) )
        {
            mexErrMsgTxt( "Error during initialization!" );
        }

        /* rfctx.residue_final_method = 0; */

        RFC_feed( &rfctx, data, data_len, 1 /*do_finalize*/  );

        if( plhs )
        {
            plhs[0] = mxCreateDoubleScalar( rfctx.pseudo_damage );

            if( nlhs > 1 && rfctx.residue )
            {
                mxArray* re = mxCreateDoubleMatrix( rfctx.residue_cnt, 1, mxREAL );
                if( re )
                {
                    size_t i;
                    double *val = mxGetPr(re);

                    for( i = 0; i < rfctx.residue_cnt; i++ )
                    {
                        *val++ = (double)rfctx.residue[i].value;
                    }
                    plhs[1] = re;
                }
            }

            if( nlhs > 2 && rfctx.matrix )
            {
                mxArray* matrix = mxCreateDoubleMatrix( class_count, class_count, mxREAL );
                if( matrix )
                {
                    memcpy( mxGetPr(matrix), rfctx.matrix, sizeof(double) * class_count * class_count );
                    plhs[2] = matrix;
                }
            }
            
            if( nlhs > 3 && rfctx.rp )
            {
                mxArray* rp = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( rp )
                {
                    memcpy( mxGetPr(rp), rfctx.rp, sizeof(double) * class_count );
                    plhs[3] = rp;
                }
            }

            if( nlhs > 4 && rfctx.lc )
            {
                mxArray* lc = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( lc )
                {
                    memcpy( mxGetPr(lc), rfctx.lc, sizeof(double) * class_count );
                    plhs[4] = lc;
                }
            }

            if( nlhs > 5 && rfctx.tp )
            {
                mxArray* tp = mxCreateDoubleMatrix( rfctx.tp_cnt, 2, mxREAL );
                if( tp )
                {
                    size_t i;
                    double *idx = mxGetPr(tp) + 0;
                    double *val = mxGetPr(tp) + rfctx.tp_cnt;

                    for( i = 0; i < rfctx.tp_cnt; i++ )
                    {
                        *val++ = (double)rfctx.tp[i].value;
                        *idx++ = (double)rfctx.tp[i].pos;
                    }
                    plhs[5] = tp;
                }
            }

        }

        RFC_deinit( &rfctx );
    }
}
#endif
