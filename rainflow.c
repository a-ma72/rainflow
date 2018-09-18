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
"[pd,rm] = rfc( data, class_count, class_width, class_offset, hysteresis )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"
#pragma message(RFC_MEX_USAGE)
#include <string.h>
#include <mex.h>
#endif

/* Core */
static void             RFC_feed_handle_tp          ( rfctx_s *rfctx, value_tuple_s* tp );
static void             RFC_feed_finalize           ( rfctx_s* rfctx );
static value_tuple_s *  RFC_tp_next_default         ( rfctx_s *, value_tuple_s *pt );
static void             RFC_cycle_find_4ptm         ( rfctx_s * );
static void             RFC_cycle_process           ( rfctx_s *, value_tuple_s *from, value_tuple_s *to, int flags );
/* Residual methods */
static bool             RFC_finalize_default        ( rfctx_s * );
/* Other */
static double           RFC_damage_calc_default     ( rfctx_s *, unsigned class_from, unsigned class_to );
static RFC_value_type   value_delta                 ( RFC_value_type from, RFC_value_type to, int *sign_ptr );


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
 *
 * @return     true on success
 */
bool RFC_init( void *ctx, unsigned class_count, RFC_value_type class_width, RFC_value_type class_offset, 
                          RFC_value_type hysteresis )
{
    rfctx_s       *rfctx = (rfctx_s*)ctx;
    value_tuple_s  nil   = { 0.0 };  /* All other members are zero-initialized, see ISO/IEC 9899:TC3, 6.7.8 (21) */

    if( !rfctx || rfctx->state != RFC_STATE_INIT0 ) return false;

    assert( rfctx->version == sizeof(rfctx_s) );

    /* Flags */
    rfctx->flags = RFC_FLAGS_COUNT_ALL;
    
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

    /* Memory allocation functions */
    rfctx->mem_alloc            = calloc;
    rfctx->mem_free             = free;

    /* Residue */
    rfctx->residue_cnt          = 0;
    rfctx->residue_size         = 2 * rfctx->class_count; /* max size is 2*n-1 plus interim point = 2*n */
    rfctx->residue              = (value_tuple_s*)rfctx->mem_alloc( rfctx->residue_size, sizeof(value_tuple_s) );

    /* Non-sparse storages (optional, may be NULL) */
    rfctx->matrix               = (RFC_counts_type*)rfctx->mem_alloc( class_count * class_count, sizeof(RFC_value_type) );

    /* Damage */
    rfctx->pseudo_damage        = 0.0;

    rfctx->internal.slope       = 0;
    rfctx->internal.extrema[0]  = nil;
    rfctx->internal.extrema[1]  = nil;

    if( !rfctx->residue || !rfctx->matrix )
    {
        RFC_deinit( rfctx );
        return false;
    }
    
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

    rfctx->residue             = NULL;
    rfctx->residue_size        = 0;
    rfctx->residue_cnt         = 0;

    rfctx->matrix              = NULL;

    rfctx->internal.slope      = 0;
    rfctx->internal.extrema[0] = nil;
    rfctx->internal.extrema[1] = nil;
    rfctx->internal.pos        = 0;

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

        RFC_feed_handle_tp( rfctx, &tp );
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
 */
static
void RFC_feed_handle_tp( rfctx_s *rfctx, value_tuple_s* tp )
{
    /* Test if a new turning point exists */
    if( RFC_tp_next_default( rfctx, tp ) )
    {
        /* New turning point, check for closed cycles and count */
        RFC_cycle_find_4ptm( rfctx );
    }
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
    rfctx->state = RFC_finalize_default( rfctx ) ? RFC_STATE_FINISHED : RFC_STATE_ERROR;
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
        /* Check if a new cycle is closed now */
        RFC_cycle_find_4ptm( rfctx );
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

        D_i = exp( fabs(k)  * ( log(Sa_i) - SD_log ) - ND_log );
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
 * @brief      Add one data sample (turning point) to the residue.
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

    slope  = rfctx->internal.slope;

    if( !rfctx->residue_cnt )
    {
        /* Residue is empty */

        int i_first_point;  /* Index to internal.extrema indicating 1st point */

        if( rfctx->state == RFC_STATE_INIT )
        {
            /* Very first point, initialize local min-max search */
            rfctx->internal.extrema[0] = rfctx->internal.extrema[1] = *pt;
            rfctx->state               = RFC_STATE_BUSY;
        }
        else
        {
            /* Still searching for first turning point */

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

        /* Hysteresis Filtering */
        delta = value_delta( rfctx->internal.extrema[0].value, rfctx->internal.extrema[1].value, NULL /* sign_ptr */ );

        if( delta > rfctx->hysteresis )
        {
            /* Criteria met, first turning point found */
            assert( rfctx->residue_cnt < rfctx->residue_size );
            rfctx->residue[rfctx->residue_cnt++] = rfctx->internal.extrema[i_first_point];

            /* Append interim turning point */
            assert( rfctx->residue_cnt < rfctx->residue_size );
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
        assert( rfctx->residue_cnt < rfctx->residue_size );
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

    /* Quantize "from" */
    class_from = QUANTIZE( rfctx, from->value );

    if( class_from >= rfctx->class_count ) class_from = rfctx->class_count - 1;

    /* Quantize "to" */
    class_to = QUANTIZE( rfctx, to->value );

    if( class_to >= rfctx->class_count ) class_to = rfctx->class_count - 1;
    
    /* class_from and class_to are base 0 now */

    /* Do several countings */
    if( class_from != class_to )
    {
        /* Cumulate pseudo damage */
        double damage;
        damage = RFC_damage_calc_default( rfctx, class_from, class_to );
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

        if( !RFC_init( &rfctx, class_count, class_width, class_offset, hysteresis ) )
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
        }

        RFC_deinit( &rfctx );
    }
}
#endif
