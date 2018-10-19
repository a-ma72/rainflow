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
#define GREATEST_FPRINTF greatest_fprintf
#include "greatest/greatest.h"
#include <stdarg.h>
#include <string.h>
#include <mex.h>
#endif

/* Core */
static void                 RFC_feed_handle_tp                  ( rfc_ctx_s *, const rfc_value_tuple_s* tp, bool is_right_margin );
static rfc_value_tuple_s *  RFC_tp_next                         ( rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 RFC_cycle_find_4ptm                 ( rfc_ctx_s * );
static void                 RFC_cycle_process                   ( rfc_ctx_s *, const rfc_value_tuple_s *from, const rfc_value_tuple_s *to, int flags );
/* Residual methods */
static void                 RFC_finalize_interim                ( rfc_ctx_s * );
static bool                 RFC_finalize_res                    ( rfc_ctx_s * );
static bool                 RFC_finalize_res_ignore             ( rfc_ctx_s * );
static bool                 RFC_finalize_res_halfcycles         ( rfc_ctx_s * );
static bool                 RFC_finalize_res_fullcycles         ( rfc_ctx_s * );
static bool                 RFC_finalize_res_clormann_seeger    ( rfc_ctx_s * );
static bool                 RFC_finalize_res_din                ( rfc_ctx_s * );
static bool                 RFC_finalize_res_repeated           ( rfc_ctx_s * );
/* Memory allocator */
static void *               RFC_mem_alloc_default               ( void *ptr, size_t num, size_t size );
/* Other */
static void                 RFC_tp_add                          ( rfc_ctx_s *, const rfc_value_tuple_s *pt );
static void                 RFC_tp_lock                         ( rfc_ctx_s *, bool do_lock );
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
 * @param[in]  tp            Pointer to turning points buffer
 * @param[in]  tp_cap        Number of turning points in buffer
 *
 * @return     true on success
 */
bool RFC_init( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                          RFC_value_type hysteresis, int residual_method,
                          rfc_value_tuple_s *tp, size_t tp_cap )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil   = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx || rfc_ctx->state != RFC_STATE_INIT0 ) return false;

    assert( rfc_ctx->version == sizeof(rfc_ctx_s) );

    /* Flags */
    rfc_ctx->flags = RFC_FLAGS_COUNT_ALL;
    
    /* Residual method */
    rfc_ctx->residual_method         = residual_method;

    /* Counter increments */
    rfc_ctx->full_inc                = RFC_FULL_CYCLE_INCREMENT;
    rfc_ctx->half_inc                = RFC_HALF_CYCLE_INCREMENT;
    rfc_ctx->curr_inc                = RFC_FULL_CYCLE_INCREMENT;

    if( !class_count || class_count > 512 ) return false;
    if( class_width <= 0.0 )                return false;

    /* Rainflow class parameters */
    rfc_ctx->class_count             = class_count;
    rfc_ctx->class_width             = class_width;
    rfc_ctx->class_offset            = class_offset;
    rfc_ctx->hysteresis              = ( hysteresis < 0.0 ) ? class_width : hysteresis;

    /* Woehler curve (fictive) */
    rfc_ctx->wl_sd                   =  1e3;            /* Fictive amplitude */
    rfc_ctx->wl_nd                   =  1e7;            /* Fictive count */
    rfc_ctx->wl_k                    = -5.0;            /* Fictive gradient */
    rfc_ctx->wl_k2                   =  rfc_ctx->wl_k;  /* "Miner elementar", if k == k2 */
    rfc_ctx->wl_omission             =  0.0;            /* No omission per default */

    /* Memory allocator */
    if( !rfc_ctx->mem_alloc )
    {
      rfc_ctx->mem_alloc = RFC_mem_alloc_default;
    }
    
#if RFC_USE_DELEGATES
    /* Delegates (optional, set to NULL for standard or only your own functions! ) */
    rfc_ctx->tp_next_fcn             = NULL;
    rfc_ctx->tp_add_fcn              = NULL;
    rfc_ctx->cycle_find_fcn          = NULL;
    rfc_ctx->finalize_fcn            = NULL;
    rfc_ctx->damage_calc_fcn         = NULL;
#endif

    /* Residue */
    rfc_ctx->residue_cnt             = 0;
    rfc_ctx->residue_cap             = 2 * rfc_ctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */
    rfc_ctx->residue                 = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, rfc_ctx->residue_cap, sizeof(rfc_value_tuple_s) );

    /* Non-sparse storages (optional, may be NULL) */
    rfc_ctx->matrix                  = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count * class_count, sizeof(RFC_value_type) );
    rfc_ctx->rp                      = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count,               sizeof(RFC_value_type) );
    rfc_ctx->lc                      = (RFC_counts_type*)rfc_ctx->mem_alloc( NULL, class_count,               sizeof(RFC_value_type) );

    /* Damage */
    rfc_ctx->pseudo_damage           = 0.0;

    rfc_ctx->internal.slope          = 0;
    rfc_ctx->internal.extrema[0]     = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]     = nil;  /* local maximum */
    rfc_ctx->internal.margin[0]      = nil;  /* left  margin */
    rfc_ctx->internal.margin[1]      = nil;  /* right margin */
    rfc_ctx->internal.tp_delayed     = nil;

    if( !rfc_ctx->residue || !rfc_ctx->matrix || !rfc_ctx->rp || !rfc_ctx->lc )
    {
        RFC_deinit( rfc_ctx );
        return false;
    }
    
    /* Turning points storage (optional, may be NULL) */
    rfc_ctx->tp                   = tp;
    rfc_ctx->tp_cap               = tp ? tp_cap : 0;
    rfc_ctx->tp_cnt               = 0;
    rfc_ctx->tp_locked            = 0;

    rfc_ctx->state = RFC_STATE_INIT;
    return true;
}


/**
 * @brief      De-initialization (freeing memory).
 *             Turning points must be freed before!
 *
 * @param      ctx  The rainflow context
 */
void RFC_deinit( void *ctx )
{
    rfc_ctx_s         *rfc_ctx = (rfc_ctx_s*)ctx;
    rfc_value_tuple_s  nil     = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfc_ctx ) return;

    assert( rfc_ctx->version == sizeof( rfc_ctx_s ) );

    if( rfc_ctx->residue )           rfc_ctx->mem_alloc( rfc_ctx->residue, 0, 0 );
    if( rfc_ctx->matrix )            rfc_ctx->mem_alloc( rfc_ctx->matrix,  0, 0 );
    if( rfc_ctx->rp )                rfc_ctx->mem_alloc( rfc_ctx->rp,      0, 0 );
    if( rfc_ctx->lc )                rfc_ctx->mem_alloc( rfc_ctx->lc,      0, 0 );
    if( rfc_ctx->tp )                rfc_ctx->mem_alloc( rfc_ctx->tp,      0, 0 );

    rfc_ctx->residue                 = NULL;
    rfc_ctx->residue_cap             = 0;
    rfc_ctx->residue_cnt             = 0;

    rfc_ctx->matrix                  = NULL;
    rfc_ctx->rp                      = NULL;
    rfc_ctx->lc                      = NULL;

    rfc_ctx->internal.slope          = 0;
    rfc_ctx->internal.extrema[0]     = nil;  /* local minimum */
    rfc_ctx->internal.extrema[1]     = nil;  /* local maximum */
    rfc_ctx->internal.margin[0]      = nil;  /* left margin */
    rfc_ctx->internal.margin[1]      = nil;  /* right margin */
    rfc_ctx->internal.pos            = 0;
    rfc_ctx->internal.tp_delayed     = nil;

    rfc_ctx->tp                      = NULL;
    rfc_ctx->tp_cap                  = 0;
    rfc_ctx->tp_cnt                  = 0;
    rfc_ctx->tp_locked               = 0;

    rfc_ctx->state = RFC_STATE_INIT0;
}


/**
 * @brief      "Feed" counting algorithm with data samples.
 *
 * @param      ctx          The rainflow context
 * @param[in]  data         The data
 * @param[in]  data_count   The data count
 */
void RFC_feed( void *ctx, const RFC_value_type * data, size_t data_count )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx ) return;

    assert( !data_count || data );
    assert( rfc_ctx->version == sizeof( rfc_ctx_s ) );

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        assert( false );
        return;
    }

    /* Process data */
    while( data_count-- )
    {
        rfc_value_tuple_s tp = { (RFC_value_type)*data++ };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

        /* Assign class and global position (base 1) */
        tp.class = QUANTIZE( rfc_ctx, tp.value );
        tp.pos   = ++rfc_ctx->internal.pos;

        RFC_feed_handle_tp( rfc_ctx, &tp, !data_count /* is_right_margin */ );
    }
}


/**
 * @brief      Feed counting algorithm with data tuples.
 *
 * @param      ctx          The rainflow context
 * @param[in]  data         The data tuples
 * @param[in]  data_count   The data count
 */
void RFC_feed_tuple( void *ctx, const rfc_value_tuple_s *data, size_t data_count )
{
    rfc_ctx_s *rfc_ctx = (rfc_ctx_s*)ctx;

    if( !rfc_ctx ) return;

    assert( !data_count || data );
    assert( rfc_ctx->version == sizeof( rfc_ctx_s ) );

    if( rfc_ctx->state < RFC_STATE_INIT || rfc_ctx->state >= RFC_STATE_FINISHED )
    {
        assert( false );
        return;
    }

    /* Process data */
    while( data_count-- )
    {
        RFC_feed_handle_tp( rfc_ctx, data++, !data_count /* is_right_margin */ );
    }
}


/**
 * @brief      Processing data (one data point).
 *             Find turning points and check for closed cycles.
 *             Locks turning points queue, if is_right_margin.
 *
 * @param      rfc_ctx          The rainflow context
 * @param[in]  pt               The data tuple
 * @param[in]  is_right_margin  True, if this is the last data tuple
 */
static
void RFC_feed_handle_tp( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s* pt, bool is_right_margin )
{
    rfc_value_tuple_s *tp_residue;
    rfc_value_tuple_s  tp_delayed;
    bool               do_margin;

    /* Check for next turning point */
    tp_residue = RFC_tp_next( rfc_ctx, pt );

    /* Delay stage when RFC_FLAGS_ENFORCE_MARGIN is set */
    do_margin = (rfc_ctx->flags & RFC_FLAGS_ENFORCE_MARGIN) && !rfc_ctx->tp_locked;
    if( do_margin )
    {
        /* Check for left margin */
        if( rfc_ctx->internal.pos == 1 )
        {
            /* Save left margin */
            rfc_ctx->internal.margin[1]  =
            rfc_ctx->internal.tp_delayed = *pt;
            tp_residue = NULL;
        }
        else if( tp_residue && rfc_ctx->internal.tp_delayed.value == tp_residue->value )
        {
            if( rfc_ctx->internal.tp_delayed.pos > 1 )
            {
                rfc_ctx->internal.tp_delayed = *tp_residue;
                tp_residue = NULL;
            }
            else
            {
                tp_residue = NULL;
            }
        }

        if( tp_residue )
        {
            /* Emit delayed turning point */
            tp_delayed = rfc_ctx->internal.tp_delayed;
            rfc_ctx->internal.tp_delayed = *tp_residue;
            tp_residue = &tp_delayed;
        }
    }

    /* Add turning point and check for closed cycles */
    for( int i = 0; i < 2; i++ )
    {
        if( tp_residue )
        {
            /* Add new turning point */
            RFC_tp_add( rfc_ctx, tp_residue );

            /* New turning point, check for closed cycles and count */
            RFC_cycle_find_4ptm( rfc_ctx );
        }

        if( is_right_margin && do_margin )
        {
            /* Empty delay stage */
            tp_residue = &rfc_ctx->internal.tp_delayed;
        }
        else break;
    }

    if( is_right_margin )
    {
        /* Save right margin */
        rfc_ctx->internal.margin[1] = *pt;

        if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
        {
            rfc_value_tuple_s *tp_interim = &rfc_ctx->residue[rfc_ctx->residue_cnt];

            if( !do_margin || tp_interim->value != pt->value )
            {
                RFC_tp_add( rfc_ctx, tp_interim );
            }
        }
 
        if( do_margin )
        {
            RFC_tp_add( rfc_ctx, pt );
        }
    }
}


/**
 * @brief      Finalizing.
 *             Take residue into account.
 *
 * @param      rfc_ctx        The rainflow context
 */
void RFC_feed_finalize( rfc_ctx_s* rfc_ctx )
{
    /* Finalizing (process last turning point */
    rfc_ctx->state = RFC_finalize_res( rfc_ctx ) ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
}


/**
 * @brief      Handling interim turning point.
 *             There is still one unhandled turning point left.
 *             "Finalizing" takes this into account for the rainflow algorithm.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
void RFC_finalize_interim( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    /* Adjust residue: Incorporate interim turning point */
    if( rfc_ctx->state == RFC_STATE_BUSY_INTERIM )
    {
        assert( rfc_ctx->residue && rfc_ctx->residue_cnt );

        rfc_ctx->residue_cnt++;

        /* Check if a new cycle is closed now */
        RFC_cycle_find_4ptm( rfc_ctx );
    }

    rfc_ctx->state = RFC_STATE_FINALIZE;
}


/**
 * @brief      Finalize pending counts, ignoring residue.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

#if RFC_USE_DELEGATES
    if( rfc_ctx->finalize_fcn )
    {
        return rfc_ctx->finalize_fcn( rfc_ctx );
    }
#endif

    switch( rfc_ctx->residual_method )
    {
        case RFC_RES_NONE:
        case RFC_RES_IGNORE:
            return RFC_finalize_res_ignore( rfc_ctx );
        case RFC_RES_HALFCYCLES:
            return RFC_finalize_res_halfcycles( rfc_ctx );
        case RFC_RES_FULLCYCLES:
            return RFC_finalize_res_fullcycles( rfc_ctx );
        case RFC_RES_CLORMANN_SEEGER:
            return RFC_finalize_res_clormann_seeger( rfc_ctx );
        case RFC_RES_REPEATED:
            return RFC_finalize_res_repeated( rfc_ctx );
        case RFC_RES_RP_DIN:
            return RFC_finalize_res_din( rfc_ctx );
        default:
            assert( false );
            return false;
    }
}


/**
 * @brief      Finalize pending counts, ignoring residue.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_ignore( rfc_ctx_s *rfc_ctx )
{
    /* Just include interim turning point */
    return RFC_finalize_interim( rfc_ctx );
}


/**
 * @brief      Finalize pending counts, half cycles method.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_halfcycles( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
    {
    }

    return true;
}


/**
 * @brief      Finalize pending counts, full cycles method.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_fullcycles( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
    {
    }

    return true;
}


/**
 * @brief      Finalize pending counts, Clormann/Seeger method.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_clormann_seeger( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
    {
    }

    return true;
}


/**
 * @brief      Finalize pending counts, DIN method.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_din( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
    {
    }

    return true;
}


/**
 * @brief      Finalize pending counts, repeated residue method.
 *
 * @param      rfc_ctx  The rainflow context
 */
static
bool RFC_finalize_res_repeated( rfc_ctx_s *rfc_ctx )
{
    assert( rfc_ctx && rfc_ctx->state < RFC_STATE_FINALIZE );

    if( rfc_ctx->residue && rfc_ctx->residue_cnt )
    {
        /* Include interim turning point as new data series, 
           but don't modify residue history itself */
        size_t              cnt     = rfc_ctx->residue_cnt;
        rfc_value_tuple_s  *residue = (rfc_value_tuple_s*)rfc_ctx->mem_alloc( NULL, ++cnt, sizeof(rfc_value_tuple_s) );

        if( residue )
        {
            /* Make a copy of the residue */
            size_t n = cnt;
            const rfc_value_tuple_s *from = rfc_ctx->residue;
                  rfc_value_tuple_s *to   = residue;

            while( n-- )
            {
                *to++ = *from++;
            }

            /* Make last turning point interim again */
            rfc_ctx->residue_cnt--;
            
            /* Feed again with the copy */
            RFC_feed_tuple( rfc_ctx, residue, cnt );

            /* Include interim turning point again */
            rfc_ctx->residue_cnt++;

            rfc_ctx->mem_alloc( residue, 0, 0 );
        }
        else return false;

        /* Handle interim turning point */
        return RFC_finalize_interim( rfc_ctx );
    }
    return true;
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

#if RFC_USE_DELEGATES
    if( rfc_ctx->damage_calc_fcn )
    {
        return rfc_ctx->damage_calc_fcn( rfc_ctx, class_from, class_to );
    }
#endif

    /* Constants for the woehler curve */
    const double SD_log  = log(rfc_ctx->wl_sd);
    const double ND_log  = log(rfc_ctx->wl_nd);
    const double k       = rfc_ctx->wl_k;
    const double k2      = rfc_ctx->wl_k2;
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

        if( Sa_i > rfc_ctx->wl_omission )
        {
            if( Sa_i > rfc_ctx->wl_sd )
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
 * @param[in]   from      Left hand value
 * @param[in]   to        Right hand value
 * @param[out]  sign_ptr  Pointer to catch sign (may be NULL)
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

    assert( rfc_ctx && pt );
    assert( rfc_ctx->state >= RFC_STATE_INIT && rfc_ctx->state <= RFC_STATE_BUSY_INTERIM );

#if RFC_USE_DELEGATES
    if( rfc_ctx->tp_next_fcn )
    {
        return rfc_ctx->tp_next_fcn( rfc_ctx, pt );
    }
#endif

    slope = rfc_ctx->internal.slope;

    /* Handle first turning point(s) */
    if( rfc_ctx->state < RFC_STATE_BUSY_INTERIM )
    {
        /* Residue is empty, still searching first turning point(s) */

        if( rfc_ctx->state == RFC_STATE_INIT )
        {
            /* Very first point, initialize local min-max search, save left margin */
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

                slope = is_falling_slope ? -1 : 1;

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
        }
        else if( pt->value > rfc_ctx->internal.extrema[1].value )
        {
            /* Maximum */
            rfc_ctx->internal.extrema[1] = *pt;
        }
#endif

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
            do_append = 1;
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

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->cycle_find_fcn )
    {
        rfc_ctx->cycle_find_fcn( rfc_ctx );
    }
    else
#endif
    {
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

    /* If flag RFC_FLAGS_ENFORCE_MARGIN is set, cycles less than hysteresis are possible */
    if( flags & RFC_FLAGS_ENFORCE_MARGIN )
    {
        if( value_delta( from->value, to->value, NULL /* sign_ptr */ ) <= rfc_ctx->hysteresis )
        {
            return;
        }
    }

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
    }
}


/**
 * @brief   Append one turning point to the queue
 *
 * @param       rfc_ctx     The rainflow context
 * @param[in]   tp          New turning points
 */
static
void RFC_tp_add( rfc_ctx_s *rfc_ctx, const rfc_value_tuple_s *tp )
{
    assert( rfc_ctx );

    /* Add new turning point */

#if RFC_USE_DELEGATES
    /* Check for delegates */
    if( rfc_ctx->tp_add_fcn )
    {
        /* Add turning point */
        rfc_ctx->tp_add_fcn( rfc_ctx, tp, do_lock );
    }
    else
#endif
    {
        if( rfc_ctx->tp && !rfc_ctx->tp_locked )
        {
            if( tp )
            {
                assert( rfc_ctx->tp_cnt < rfc_ctx->tp_cap );
                rfc_ctx->tp[ rfc_ctx->tp_cnt++ ] = *tp;
            }
        }
    }
}


/**
 * @brief   Lock turning points queue
 *
 * @param       rfc_ctx     The rainflow context
 * @param[in]   do_lock     Turning point storage will be locked, if true
 */
static 
void RFC_tp_lock( rfc_ctx_s *rfc_ctx, bool do_lock )
{
    rfc_ctx->tp_locked = do_lock;
}


/**
 * @brief       (Re-)Allocate or free memory
 * @param       ptr     Previous data pointer, or NULL, if unset
 * @param[in]   num     The number of elements
 * @param[in]   size    The size of one element in bytes
 * @returns     New memory pointer or NULL if either num or size is 0
 * 
 */
static
void * RFC_mem_alloc_default( void *ptr, size_t num, size_t size )
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


void RFC_lc_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_len );
void RFC_rp_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_len );

/**
 * Calculate level crossing counts from rainflow matrix, write results to buffer.
 */
void RFC_lc_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_len )
{
    unsigned i, j, k;
    bool     up     = rfc_ctx->flags & RFC_FLAGS_COUNT_LC_UP;
    bool     dn     = rfc_ctx->flags & RFC_FLAGS_COUNT_LC_DN;
    size_t   maxcnt = buffer_len / sizeof(RFC_counts_type);

    assert( rfc_ctx );

    if( !buffer || !maxcnt ) return;

    for( i = 0; i < rfc_ctx->class_count; i++ ) 
    {
        RFC_counts_type counts = (RFC_counts_type)0;

        for( j = i; j < rfc_ctx->class_count; j++ )  /* To */
        {
            for( k = 0; k < i; k++ )               /* From */
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

        buffer[i > maxcnt ? maxcnt : i] = counts;
    }
}


/**
 * Calculate range pair counts from rainflow matrix, write results to buffer.
 */
void RFC_rp_from_matrix( rfc_ctx_s *rfc_ctx, RFC_counts_type* buffer, size_t buffer_len )
{
    unsigned i, j;
    size_t   maxcnt = buffer_len / sizeof(RFC_counts_type);

    assert( rfc_ctx && buffer && buffer_len >= rfc_ctx->class_count * sizeof(RFC_counts_type) );

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
}


/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


#ifdef MATLAB_MEX_FILE


int greatest_fprintf( FILE* f, const char* fmt, ... )
{
    va_list al;
    char *buffer = NULL;
    int len;

    va_start( al, fmt ); 
    len = vsnprintf( buffer, 0, fmt, al );
    if( len > 0 )
    {
        buffer = (char*)calloc( len + 1, 1 );
        if( buffer )
        {
            va_start( al, fmt );
            (void)vsnprintf( buffer, len + 1, fmt, al );
            mexPrintf( "%s", buffer );
            free( buffer );
        }
    }
    va_end( al );

    return len;
}

#define INIT_ARRAY(...) __VA_ARGS__
#define SIMPLE_RFC(TP,TP_N,OFFS,X,X_N) \
    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */,                                                  \
                        RFC_RES_NONE /* residual_method */,                                  \
                        TP /* *tp */, TP_N /* tp_cap */ ) )                                  \
    {                                                                                        \
        RFC_VALUE_TYPE data[] = {INIT_ARRAY X};                                             \
        RFC_feed( &ctx, data, X_N );                                                         \
        RFC_feed_finalize( &ctx );                                                           \
    }                                                                                        \
    else FAIL();

#define SIMPLE_RFC2(TP,TP_N,OFFS,X,X_N) \
    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */,                                                  \
                        RFC_RES_NONE /* residual_method */,                                  \
                        TP /* *tp */, TP_N /* tp_cap */ ) )                                  \
    {                                                                                        \
        RFC_VALUE_TYPE data[] = {INIT_ARRAY X};                                             \
        ctx.flags |= RFC_FLAGS_ENFORCE_MARGIN;                                               \
        RFC_feed( &ctx, data, X_N );                                                         \
        RFC_feed_finalize( &ctx );                                                           \
    }                                                                                        \
    else FAIL();


/*

*/
TEST RFC_test_turning_points(void)
{
    rfc_ctx_s         ctx = {sizeof(ctx)};
    rfc_value_tuple_s tp[10];

    SIMPLE_RFC( tp, 10, 0.0, (0), 0 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC( tp, 10, 0.0, (0), 1 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC( tp, 10, 0.0, (0,0), 2 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC( tp, 10, 0.0, (0,0.1), 2 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC( tp, 10, 0.0, (0,1), 2 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC( tp, 10, 0.0, (1.0, 1.1, 1.2, 1.1, 1.3, 1.0, 1.98, 1.0), 8 );
    ASSERT( ctx.tp_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

    SIMPLE_RFC2( tp, 10, 0.0, (1.0, 1.1, 1.2, 1.1, 1.3, 1.0, 1.98, 1.0), 8 );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( ctx.tp[0].value == 1.0 && ctx.tp[0].pos == 1 );
    ASSERT( ctx.tp[1].value == 1.0 && ctx.tp[1].pos == 8 );
    ASSERT( ctx.residue_cnt == 0 );
    ctx.tp = NULL;
    RFC_deinit( &ctx );

#if 0
    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, 0.0 /* class_offset */, 
                        1 /* hysteresis */, 
                        RFC_RES_NONE /* residual_method */,
                        tp /* *tp */, 10 /* tp_cap */ ) )
    {
        double data[] = { 1.0, 1.1, 1.2, 2.0, 2.1, 1.1, 1.3, 1.0, 1.98, 1.0 };
        RFC_feed( &ctx, data, NUMEL(data) );
        RFC_feed_finalize( &ctx );
        ASSERT( ctx.tp_cnt == 3 );
        ASSERT( ctx.tp[0].value == 1.0 && ctx.tp[0].pos == 1 );
        ASSERT( ctx.tp[1].value == 2.1 && ctx.tp[1].pos == 5 );
        ASSERT( ctx.tp[2].value == 1.0 && ctx.tp[2].pos == 8 );
        ASSERT( ctx.residue_cnt == 3 );
        ASSERT( ctx.residue[0].value == 1.0 && ctx.residue[0].pos == 1 );
        ASSERT( ctx.residue[1].value == 2.1 && ctx.residue[1].pos == 5 );
        ASSERT( ctx.residue[2].value == 1.0 && ctx.residue[2].pos == 8 );
        ctx.tp = NULL;
        RFC_deinit( &ctx );
    }
    else FAIL();

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, 0.0 /* class_offset */, 
                        1 /* hysteresis */, 
                        RFC_RES_NONE /* residual_method */,
                        tp /* *tp */, 10 /* tp_cap */ ) )
    {
        double data[] = { 1.0, 2.1, 1.0, 1.0 };
        ctx.flags = RFC_FLAGS_ENFORCE_MARGIN;
        RFC_feed( &ctx, data, NUMEL(data) );
        RFC_feed_finalize( &ctx );
        ASSERT( ctx.tp_cnt == 4 );
        ASSERT( ctx.tp[0].value == 1.0 && ctx.tp[0].pos == 1 );
        ASSERT( ctx.tp[1].value == 2.1 && ctx.tp[1].pos == 5 );
        ASSERT( ctx.tp[2].value == 1.0 && ctx.tp[2].pos == 8 );
        ASSERT( ctx.tp[3].value == 1.0 && ctx.tp[3].pos == 10 );
        ASSERT( ctx.residue_cnt == 4 );
        ASSERT( ctx.residue[0].value == 1.0 && ctx.residue[0].pos == 1 );
        ASSERT( ctx.residue[1].value == 2.1 && ctx.residue[1].pos == 5 );
        ASSERT( ctx.residue[2].value == 1.0 && ctx.residue[2].pos == 8 );
        ASSERT( ctx.residue[3].value == 1.0 && ctx.residue[3].pos == 10 );
        ctx.tp = NULL;
        RFC_deinit( &ctx );
    }
    else FAIL();
#endif

    PASS();
}

/* local suite (greatest) */
SUITE(RFC_TURNING_POINTS)
{
    RUN_TEST( RFC_test_turning_points );
}

GREATEST_MAIN_DEFS();

int RFC_test_main( int argc, char* argv[] )
{
    GREATEST_MAIN_BEGIN();      /* init & parse command-line args */
    RUN_SUITE( RFC_TURNING_POINTS );
    GREATEST_MAIN_END();        /* display results */        
}

/**
 * MATLAB wrapper for the rainflow algorithm
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
    if( !nrhs )
    {
        mexPrintf( "%s", RFC_MEX_USAGE );
        mexPrintf( "%s\n", "Running self tests..." );

        RFC_test_main( 0, NULL );

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

        double   *data         = mxGetPr( mxData );
        size_t    data_len     = mxGetNumberOfElements( mxData );
        unsigned  class_count  = (unsigned)( mxGetScalar( mxClassCount ) + 0.5 );
        double    class_width  = mxGetScalar( mxClassWidth );
        double    class_offset = mxGetScalar( mxClassOffset );
        double    hysteresis   = mxGetScalar( mxHysteresis );

        rfc_ctx.tp = (rfc_value_tuple_s*)calloc( data_len, sizeof(rfc_value_tuple_s) );

        if( !RFC_init( &rfc_ctx, class_count, class_width, class_offset, hysteresis, RFC_RES_NONE, rfc_ctx.tp, data_len ) )
        {
            mexErrMsgTxt( "Error during initialization!" );
        }

        /* rfc_ctx.residue_final_method = 0; */

        RFC_feed( &rfc_ctx, data, data_len  );

        if( plhs )
        {
            plhs[0] = mxCreateDoubleScalar( rfc_ctx.pseudo_damage );

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

            if( nlhs > 2 && rfc_ctx.matrix )
            {
                mxArray* matrix = mxCreateDoubleMatrix( class_count, class_count, mxREAL );
                if( matrix )
                {
                    memcpy( mxGetPr(matrix), rfc_ctx.matrix, sizeof(double) * class_count * class_count );
                    plhs[2] = matrix;
                }
            }
            
            if( nlhs > 3 && rfc_ctx.rp )
            {
                mxArray* rp = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( rp )
                {
                    memcpy( mxGetPr(rp), rfc_ctx.rp, sizeof(double) * class_count );
                    plhs[3] = rp;
                }
            }

            if( nlhs > 4 && rfc_ctx.lc )
            {
                mxArray* lc = mxCreateDoubleMatrix( class_count, 1, mxREAL );
                if( lc )
                {
                    memcpy( mxGetPr(lc), rfc_ctx.lc, sizeof(double) * class_count );
                    plhs[4] = lc;
                }
            }

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

        }

        RFC_deinit( &rfc_ctx );
    }
}
#endif
