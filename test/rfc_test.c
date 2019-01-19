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
 *    Test suite
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

#define GREATEST_FPRINTF fprintf
#define RFC_VALUE_TYPE   double

#include "../rainflow.h"
#include "../greatest/greatest.h"
#include <locale.h>
#include <math.h>
#include <float.h>


#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define NUMEL(x) (sizeof(x)/sizeof((x)[0]))

typedef struct mem_chunk
{
    size_t             size, 
                       count;
    struct  mem_chunk *next;
    RFC_VALUE_TYPE     data[1];
} mem_chunk;

      rfc_ctx_s   ctx              = { sizeof(ctx) };
      mem_chunk  *mem_chain        = NULL;
const char       *long_series_file = NULL;


static
mem_chunk* new_chunk( size_t size )
{
    if( !size ) return NULL;

    mem_chunk* chunk = (mem_chunk*)calloc( 1, (size-1) * sizeof(RFC_VALUE_TYPE) + sizeof(mem_chunk) );
    if( chunk )
    {
        chunk->size  = size;
        chunk->count = 0;
    }

    return chunk;
}


#if RFC_TP_SUPPORT
void export_tp( const char *filename, rfc_value_tuple_s* data, size_t count )
{
    FILE* fid = fopen( filename, "wt" );

    if( fid )
    {
        while( count-- )
        {
            //fprintf(fid, "%g\t%lu\t%lu\t%lu\t%g\n", data->value, data->class, data->tp_pos, data->pos, data->damage);
            fprintf(fid, "%g\t%d\t%llu\t%llu\n", 
                          data->value, data->cls, 
                         (long long unsigned int)data->tp_pos, (long long unsigned int)data->pos );
            data++;
        }

        fclose( fid );
    }
}
#endif /*RFC_TP_SUPPORT*/


double rfm_peek( rfc_ctx_s *rfc_ctx, int from, int to )
{
    RFC_counts_type counts;

    counts = rfc_ctx->rfm[ rfc_ctx->class_count * (from-1) + (to-1) ];
    return (double)counts / rfc_ctx->full_inc;
}


#if RFC_TP_SUPPORT
#define INIT_ARRAY(...) __VA_ARGS__
#define SIMPLE_RFC_0(CCNT,TP,TP_N,OFFS) \
    if( RFC_init( &ctx, CCNT /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) &&                             \
        RFC_tp_init( &ctx, TP /* *tp */, TP_N /* tp_cap */, true /* is_static */ ) )           \
    {                                                                                          \
        RFC_VALUE_TYPE data[] = {0};                                                           \
        ASSERT( RFC_feed( &ctx, data, 0 ) );                                                   \
        ASSERT( RFC_finalize( &ctx, RFC_RES_NONE /* residual_method */ ) );                    \
    }                                                                                          \
    else FAIL();

#define INIT_ARRAY(...) __VA_ARGS__
#define SIMPLE_RFC(CCNT,TP,TP_N,OFFS,X) \
    if( RFC_init( &ctx, CCNT /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) &&                             \
        RFC_tp_init( &ctx, TP /* *tp */, TP_N /* tp_cap */, true /* is_static */ ) )           \
    {                                                                                          \
        RFC_VALUE_TYPE data[] = {INIT_ARRAY X};                                                \
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );                                         \
        ASSERT( RFC_finalize( &ctx, RFC_RES_NONE /* residual_method */ ) );                    \
    }                                                                                          \
    else FAIL();

#define SIMPLE_RFC_MARGIN_0(CCNT,TP,TP_N,OFFS) \
    if( RFC_init( &ctx, CCNT /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */, RFC_FLAGS_COUNT_ALL|RFC_FLAGS_ENFORCE_MARGIN ) &&  \
        RFC_tp_init( &ctx, TP /* *tp */, TP_N /* tp_cap */, true /* is_static */ ) )           \
    {                                                                                          \
        RFC_VALUE_TYPE data[] = {0};                                                           \
        ASSERT( RFC_feed( &ctx, data, 0 ) );                                                   \
        ASSERT( RFC_finalize( &ctx, RFC_RES_NONE /* residual_method */ ) );                    \
    }                                                                                          \
    else FAIL();

#define SIMPLE_RFC_MARGIN(CCNT,TP,TP_N,OFFS,X) \
    if( RFC_init( &ctx, CCNT /* class_count */, 1 /* class_width */, OFFS /* class_offset */,  \
                        1 /* hysteresis */, RFC_FLAGS_COUNT_ALL|RFC_FLAGS_ENFORCE_MARGIN ) &&  \
        RFC_tp_init( &ctx, TP /* *tp */, TP_N /* tp_cap */, true /* is_static */ ) )           \
    {                                                                                          \
        RFC_VALUE_TYPE data[] = {INIT_ARRAY X};                                                \
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );                                         \
        ASSERT( RFC_finalize( &ctx, RFC_RES_NONE /* residual_method */ ) );                    \
    }                                                                                          \
    else FAIL();



TEST RFC_tp_prune_test( int ccnt )
{
    RFC_VALUE_TYPE      data[10000];
    size_t              data_len            =  NUMEL( data );
    RFC_VALUE_TYPE      x_max;
    RFC_VALUE_TYPE      x_min;
    unsigned            class_count         =  ccnt ? 100 : 0;
    RFC_VALUE_TYPE      class_width;
    RFC_VALUE_TYPE      class_offset;
    RFC_VALUE_TYPE      hysteresis;
    rfc_value_tuple_s   tp[10000]           = {0};
    size_t              i;

    if(1)
    {
#include "long_series.c"

        for( i = 0; i < data_len; i++ )
        {
            double value = data_export[i];
            data[i] = value;
            if( !i )
            {
                x_max = x_min = value;
            }
            else
            {
                if( value > x_max ) x_max = value;
                if( value < x_min ) x_min = value;
            }
        }
    }

    class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    class_offset    =  x_min - class_width / 2;
    hysteresis      =  class_width;

    GREATEST_FPRINTF( GREATEST_STDOUT, "\nTest long series:" );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass count  = %d", class_count );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass width  = %g", class_width );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass offset = %g", class_offset );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\n" );

    if( class_count )
    {
        ASSERT( class_width > 0.0 );
        ASSERT( class_count > 1 );
        ASSERT( x_min >= class_offset );
        ASSERT( x_max <  class_offset + class_width * class_count );
    }

    ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
    ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
    ASSERT( RFC_feed( &ctx, data, /* count */ data_len ) );
    ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );
    ctx.tp_locked = 0;;
    ASSERT( RFC_tp_prune( &ctx, /*count*/ 100, /*flags*/ RFC_FLAGS_TPPRUNE_PRESERVE_POS | RFC_FLAGS_TPPRUNE_PRESERVE_RES ) );

    ASSERT( ctx.tp_cnt == ccnt ? 107 : 100 );
    /* Should not change anything: */
    ASSERT( RFC_tp_prune( &ctx, /*count*/ 100, /*flags*/ RFC_FLAGS_TPPRUNE_PRESERVE_POS | RFC_FLAGS_TPPRUNE_PRESERVE_RES ) );
    ASSERT( ctx.tp_cnt == ccnt ? 107 : 100 );

    ASSERT( RFC_tp_prune( &ctx, /*count*/ 0, /*flags*/ RFC_FLAGS_TPPRUNE_PRESERVE_POS | RFC_FLAGS_TPPRUNE_PRESERVE_RES ) );
    ASSERT( ctx.tp_cnt == ctx.residue_cnt );

    for( i = 0; i < ctx.residue_cnt; i++ )
    {
        ASSERT( ctx.tp[i].tp_pos == 0 );
        ASSERT( ctx.residue[i].tp_pos == i + 1 );
        ASSERT( ctx.tp[i].damage == 0.0 );
        ASSERT( ctx.residue[i].damage == 0.0 );
        ASSERT( ctx.tp[i].pos == ctx.residue[i].pos );
        ASSERT( ctx.tp[i].value == ctx.residue[i].value );
    }

    ASSERT( RFC_tp_prune( &ctx, /*count*/ 0, /*flags*/ RFC_FLAGS_TPPRUNE_PRESERVE_POS ) );
    ASSERT( ctx.tp_cnt == 0 );

    for( i = 0; i < ctx.residue_cnt; i++ )
    {
        ASSERT( ctx.residue[i].tp_pos == 0 );
    }

    if( ctx.state != RFC_STATE_INIT0 )
    {
        ASSERT( RFC_deinit( &ctx ) );
    }

    PASS();
}


TEST RFC_test_turning_points( int ccnt )
{
    //rfc_ctx_s         ctx = {sizeof(ctx)};
    rfc_value_tuple_s tp[10];
    
    if( ccnt ) ccnt = 10;

    /*******************************************/
    /*        Test 0, 1 or 2 samples           */
    /*******************************************/
    SIMPLE_RFC_0( ccnt, tp, 10, 0.0 );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC( ccnt, tp, 10, 0.0, (0) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC( ccnt, tp, 10, 0.0, (0,0) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC( ccnt, tp, 10, 0.0, (0.0f, 0.1f) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC( ccnt, tp, 10, 0.0, (0.0f, 1.0f) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    /**************** Test margin *******************/
    SIMPLE_RFC_MARGIN_0( ccnt, tp, 10, 0.0 );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (0) );
    ASSERT( ctx.tp_cnt == 1 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (0, 0) );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (0.0f, 0.1f) );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (0.0f, 1.0f) );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( RFC_deinit( &ctx ) );

    /*******************************************/
    /*           Test longer series            */
    /*******************************************/
    /* Still in hysteresis band */
    SIMPLE_RFC( ccnt, tp, 10, 0.0, (0.0f, 0.0f, 1.0f, 1.0f) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( ctx.residue_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC( ccnt, tp, 10, 0.0, (1.0f, 1.1f, 1.2f, 1.1f, 1.3f, 1.0f, 1.98f, 1.0f) );
    ASSERT( ctx.tp_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    /* Series with 3 turning points */
    SIMPLE_RFC( ccnt, tp, 10, 0.0, (1.0f, 1.1f, 1.2f, 2.0f, 2.1f, 1.1f, 1.3f, 1.0f, 1.98f, 1.0f) );
    ASSERT( ctx.tp_cnt == 3 );
    ASSERT( ctx.tp[0].value == 1.0f && ctx.tp[0].pos == 1 );
    ASSERT( ctx.tp[1].value == 2.1f && ctx.tp[1].pos == 5 );
    ASSERT( ctx.tp[2].value == 1.0f && ctx.tp[2].pos == 8 );
    if( ctx.class_count )
    {
        ASSERT( ctx.residue_cnt == 3 );
        ASSERT( ctx.residue[0].value == 1.0f && ctx.residue[0].pos == 1 && ctx.residue[0].tp_pos == 1 );
        ASSERT( ctx.residue[1].value == 2.1f && ctx.residue[1].pos == 5 && ctx.residue[1].tp_pos == 2 );
        ASSERT( ctx.residue[2].value == 1.0f && ctx.residue[2].pos == 8 && ctx.residue[2].tp_pos == 3 );
    }
    ASSERT( RFC_deinit( &ctx ) );

    /**************** Test margin *******************/
    /* Still in hysteresis band */
    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (0.0f, 0.0f, 1.0f, 1.0f) );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( ctx.tp[0].value == 0.0f && ctx.tp[0].pos == 1 );
    ASSERT( ctx.tp[1].value == 1.0f && ctx.tp[1].pos == 4 );
    ASSERT( ctx.residue_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (1.0f, 1.1f, 1.2f, 1.1f, 1.3f, 1.0f, 1.98f, 1.0f) );
    ASSERT( ctx.tp_cnt == 2 );
    ASSERT( ctx.tp[0].value == 1.0f && ctx.tp[0].pos == 1 );
    ASSERT( ctx.tp[1].value == 1.0f && ctx.tp[1].pos == 8 );
    ASSERT( ctx.residue_cnt == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    /* Series with 3 turning points */
    SIMPLE_RFC_MARGIN( ccnt, tp, 10, 0.0, (1.0f, 1.0f, 2.1f, 2.1f, 1.0f, 1.0f) );
    ASSERT( ctx.tp_cnt == 3 );
    ASSERT( ctx.tp[0].value == 1.0f && ctx.tp[0].pos == 1 );
    ASSERT( ctx.tp[1].value == 2.1f && ctx.tp[1].pos == 3 );
    ASSERT( ctx.tp[2].value == 1.0f && ctx.tp[2].pos == 6 ); /* Turning point at right margin! */
    if( ctx.class_count )
    {
        ASSERT( ctx.residue_cnt == 3 );
        ASSERT( ctx.residue[0].value == 1.0f && ctx.residue[0].pos == 1 && ctx.residue[0].tp_pos == 1 );
        ASSERT( ctx.residue[1].value == 2.1f && ctx.residue[1].pos == 3 && ctx.residue[1].tp_pos == 2 );
        ASSERT( ctx.residue[2].value == 1.0f && ctx.residue[2].pos == 5 && ctx.residue[2].tp_pos == 3 );  /* In residue, turning point at original position! */
    }
    ASSERT( RFC_deinit( &ctx ) );

    PASS();
}
#endif /*RFC_TP_SUPPORT*/


TEST RFC_empty( int ccnt )
{
    RFC_VALUE_TYPE      x_max           =  1;
    RFC_VALUE_TYPE      x_min           = -1;
    unsigned            class_count     =  ccnt ? 100 : 0;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width;
#if RFC_TP_SUPPORT
    rfc_value_tuple_s   tp[10];
#endif /*RFC_TP_SUPPORT*/
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {0};
        RFC_VALUE_TYPE sum = 0.0;

#if RFC_TP_SUPPORT
        ASSERT( NUMEL(tp) >= NUMEL(data) );
#endif /*RFC_TP_SUPPORT*/

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
#endif /*RFC_TP_SUPPORT*/
        ASSERT( RFC_feed( &ctx, data, /* count */ 0 ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i];
        }

        ASSERT_EQ( sum, 0.0 );
        ASSERT_EQ( ctx.residue_cnt, 0 );
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        ASSERT( RFC_deinit( &ctx ) );
    }

    PASS();
}


TEST RFC_cycle_up( int ccnt )
{
    RFC_VALUE_TYPE      x_max           =  4;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  ccnt ? 4 : 0;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width * 0.99;
#if RFC_TP_SUPPORT
    rfc_value_tuple_s   tp[10];
#endif /*RFC_TP_SUPPORT*/
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {1,3,2,4};
        RFC_VALUE_TYPE sum = 0.0;

#if RFC_TP_SUPPORT
        ASSERT( NUMEL(tp) >= NUMEL(data) );
#endif /*RFC_TP_SUPPORT*/

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
#endif /*RFC_TP_SUPPORT*/
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 3, 2 ), 1.0 );
            ASSERT_EQ( ctx.residue_cnt, 2 );
            ASSERT_EQ( ctx.residue[0].value, 1.0 );
            ASSERT_EQ( ctx.residue[1].value, 4.0 );
            ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
            ASSERT_EQ( ctx.residue[0].pos, 1 );
            ASSERT_EQ( ctx.residue[1].pos, 4 );
        }
#if RFC_TP_SUPPORT
        if( class_count )
        {
            ASSERT_EQ( ctx.residue[0].tp_pos, 1 );
            ASSERT_EQ( ctx.residue[1].tp_pos, 4 );
        }
        for( i = 0; i < ctx.tp_cnt; i++ )
        {
            ASSERT_EQ( ctx.tp[i].pos, i + 1 );
            ASSERT_EQ( ctx.tp[i].tp_pos, 0 );
        }
#endif /*RFC_TP_SUPPORT*/
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        ASSERT( RFC_deinit( &ctx ) );
    }

    PASS();
}


TEST RFC_cycle_down( int ccnt )
{
    RFC_VALUE_TYPE      x_max           =  4;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  ccnt ? 4 : 0;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width * 0.99;
#if RFC_TP_SUPPORT
    rfc_value_tuple_s   tp[10];
#endif /*RFC_TP_SUPPORT*/
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {4,2,3,1};
        RFC_VALUE_TYPE sum = 0.0;

#if RFC_TP_SUPPORT
        ASSERT( NUMEL(tp) >= NUMEL(data) );
#endif /*RFC_TP_SUPPORT*/

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
#endif /*RFC_TP_SUPPORT*/
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 2, 3 ), 1.0 );
            ASSERT_EQ( ctx.residue_cnt, 2 );
            ASSERT_EQ( ctx.residue[0].value, 4.0 );
            ASSERT_EQ( ctx.residue[1].value, 1.0 );
            ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
            ASSERT_EQ( ctx.residue[0].pos, 1 );
            ASSERT_EQ( ctx.residue[1].pos, 4 );
        }
#if RFC_TP_SUPPORT
        if( class_count )
        {
            ASSERT_EQ( ctx.residue[0].tp_pos, 1 );
            ASSERT_EQ( ctx.residue[1].tp_pos, 4 );
        }
        for( i = 0; i < ctx.tp_cnt; i++ )
        {
            ASSERT_EQ( ctx.tp[i].pos, i + 1 );
            ASSERT_EQ( ctx.tp[i].tp_pos, 0 );
        }
#endif /*RFC_TP_SUPPORT*/
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        ASSERT( RFC_deinit( &ctx ) );
    }

    PASS();
}


TEST RFC_small_example( int ccnt )
{
    RFC_VALUE_TYPE      x_max           =  6;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  ccnt ? (unsigned)x_max : 0;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width;
#if RFC_TP_SUPPORT
    rfc_value_tuple_s   tp[20];
#endif /*RFC_TP_SUPPORT*/
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2};
        RFC_VALUE_TYPE sum = 0.0;

#if RFC_TP_SUPPORT
        ASSERT( NUMEL(tp) >= NUMEL(data) );
#endif /*RFC_TP_SUPPORT*/

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
#endif /*RFC_TP_SUPPORT*/
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 7.0 );
            ASSERT_EQ( rfm_peek( &ctx, 5, 3 ), 2.0 );
            ASSERT_EQ( rfm_peek( &ctx, 6, 3 ), 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 1, 4 ), 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 2, 4 ), 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 1, 6 ), 2.0 );
            ASSERT_EQ( ctx.residue_cnt, 5 );
            ASSERT_EQ( ctx.residue[0].value, 2.0 );
            ASSERT_EQ( ctx.residue[1].value, 6.0 );
            ASSERT_EQ( ctx.residue[2].value, 1.0 );
            ASSERT_EQ( ctx.residue[3].value, 5.0 );
            ASSERT_EQ( ctx.residue[4].value, 2.0 );
        }
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        ASSERT( RFC_deinit( &ctx ) );
    }

    PASS();
}


TEST RFC_long_series( int ccnt )
{
    bool                need_conf           =  false;
    bool                do_result_check     =  true;
    RFC_VALUE_TYPE      data[10000];
    size_t              data_len            =  NUMEL( data );
    RFC_VALUE_TYPE      x_max;
    RFC_VALUE_TYPE      x_min;
    unsigned            class_count         =  ccnt ? 100 : 0;
    RFC_VALUE_TYPE      class_width;
    RFC_VALUE_TYPE      class_offset;
    RFC_VALUE_TYPE      hysteresis;
    rfc_value_tuple_s   tp[10000]           = {0};
    size_t              i;

    if(1)
    {
#include "long_series.c"

        for( i = 0; i < data_len; i++ )
        {
            double value = data_export[i];
            data[i] = value;
            if( !i )
            {
                x_max = x_min = value;
            }
            else
            {
                if( value > x_max ) x_max = value;
                if( value < x_min ) x_min = value;
            }
        }
    }
    else        
    {
        mem_chunk *chunk;
        const
        size_t     chunk_size = 10 * 1024;
        FILE*      file       = NULL;
        char       buf[81]    = {0};
        int        len;
        int        i;

        ASSERT( mem_chain = chunk = new_chunk( chunk_size ) );

        file = fopen( long_series_file, "rt" );
        ASSERT( file );

        data_len = 0;
        for( i = 0; !feof(file); i++ )
        {
            double value;

            if( fgets( buf, sizeof(buf), file ) )
            {
                if( !i && ( 0 == sscanf( buf, " * %n", &len ) ) && ( strlen(buf) == len ) )
                {
                    need_conf = true;
                }
                else if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) )
                {
                    if( chunk->count == chunk->size )
                    {
                        ASSERT( chunk->next = new_chunk( chunk_size ) );
                        chunk = chunk->next;
                    }

                    chunk->data[data_len%chunk_size] = value;
                    chunk->count++;
                    if( !data_len++ )
                    {
                        x_max = x_min = value;
                    }
                    else
                    {
                        if( value > x_max ) x_max = value;
                        if( value < x_min ) x_min = value;
                    }
                }
            }
        }
        fclose( file );
    }

    if( !need_conf )
    {
        class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
        class_offset    =  x_min - class_width / 2;
        hysteresis      =  class_width;
    }
    else
    {
        char buf[81];
        int len;
        double value;

        do_result_check = false;

        GREATEST_FPRINTF( GREATEST_STDOUT, "\n%s", "Test long series:" );
        GREATEST_FPRINTF( GREATEST_STDOUT, "\nMaximum found at %g", x_max );
        GREATEST_FPRINTF( GREATEST_STDOUT, "\nMinimum found at %g\n", x_min );
        GREATEST_FPRINTF( GREATEST_STDOUT, "\nEnter class parameters:" );
        GREATEST_FPRINTF( GREATEST_STDOUT, "\n" );

        class_count = 100;
        printf( "Class count (%d): ", class_count );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) && ( value > 0.0 ) )
            {
                class_count = (unsigned)( value + 0.5 );
            }
        }

        hysteresis = class_width = (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
        printf( "Class width (%g): ", class_width );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) && ( value > 0.0 ) )
            {
                hysteresis = class_width = value;
            }
        }

        class_offset  =  x_min - class_width / 2;
        printf( "Class offset (%g): ", class_offset );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) )
            {
                class_offset = value;
            }
        }
        printf( "\n" );
    }

    GREATEST_FPRINTF( GREATEST_STDOUT, "\nTest long series:" );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass count  = %d", class_count );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass width  = %g", class_width );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\nClass offset = %g", class_offset );
    GREATEST_FPRINTF( GREATEST_STDOUT, "\n" );

    if( class_count )
    {
        ASSERT( class_width > 0.0 );
        ASSERT( class_count > 1 );
        ASSERT( x_min >= class_offset );
        ASSERT( x_max <  class_offset + class_width * class_count );
    }

    ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
#if RFC_TP_SUPPORT
    ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /* is_static */ true ) );
#endif /*RFC_TP_SUPPORT*/
    
    if( mem_chain )
    {
        mem_chunk *it = mem_chain;
        while( it )
        {
            mem_chunk *next = it->next;
            ASSERT( RFC_feed( &ctx, it->data, /* count */ it->count ) );
            free( it );
            it = mem_chain = next;
        }
    }
    else
    {
        ASSERT( RFC_feed( &ctx, data, /* count */ data_len ) );
    }
    ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

    if(1)
    {
        FILE*   file = NULL;
        int     from, to, i;

        setlocale( LC_ALL, "" );

        file = fopen( "long_series_results.txt", "wt" );
        ASSERT( file );
        fprintf( file, "Class count: %d\n", (int)ctx.class_count );
        fprintf( file, "Class width:  %.5f\n", ctx.class_width );
        fprintf( file, "Class offset:  %.5f\n", ctx.class_offset );
        fprintf( file, "Damage: %g\n", ctx.damage);
        fprintf( file, "\nfrom (int base 0);to (int base 0);from (Klassenmitte);to (Klassenmitte);counts\n" );

        for( from = 0; from < (int)ctx.class_count; from++ )
        {
            for( to = 0; to < (int)ctx.class_count; to++ )
            {
                double value = (double)ctx.rfm[from * (int)ctx.class_count + to] / ctx.full_inc;

                if( value > 0.0 )
                {
                    fprintf( file, "%d;",  from );
                    fprintf( file, "%d;",  to );
                    fprintf( file, "%g;",  ctx.class_width * (0.5 + from) + class_offset );
                    fprintf( file, "%g;",  ctx.class_width * (0.5 + to  ) + class_offset );
                    fprintf( file, "%g\n", value );
                }
            }
        }
        fprintf( file, "\n\nResidue (classes base 0):\n" );
        for( i = 0; i < ctx.residue_cnt; i++ )
        {
            fprintf( file, "%s%d", i ? ", " : "", ctx.residue[i].cls );
        }

        fprintf( file, "\n" );
        fclose( file );
    }
    
    if( do_result_check && class_count )
    {
        do
        {
            RFC_VALUE_TYPE sum = 0.0;
            double damage = 0.0;

            for( i = 0; i < class_count * class_count; i++ )
            {
                sum += ctx.rfm[i] / ctx.full_inc;
            }

            /* Check matrix sum */
            ASSERT_EQ( sum, 602.0 );
            /* Check damage value */
            GREATEST_ASSERT_IN_RANGE( ctx.damage, 4.8703e-16, 0.00005e-16 );
#if !RFC_MINIMAL
            /* Damage must equal to damage calculated in postprocess */
            ASSERT( RFC_damage_from_rp( &ctx, NULL /*rp*/, NULL /*Sa*/, &damage, RFC_RP_DAMAGE_CALC_TYPE_DEFAULT ) );
            GREATEST_ASSERT_IN_RANGE( damage, 4.8703e-16, 0.00005e-16 );
            damage = 0.0;
            ASSERT( RFC_damage_from_rfm( &ctx, NULL /*rfm*/, &damage ) );
            GREATEST_ASSERT_IN_RANGE( damage, 4.8703e-16, 0.00005e-16 );
#endif /*!RFC_MINIMAL*/
            /* Check residue */
            ASSERT_EQ( ctx.residue_cnt, 10 );
            ASSERT_EQ_FMT( ctx.residue[0].value,   0.54, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[1].value,   2.37, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[2].value,  -0.45, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[3].value,  17.45, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[4].value, -50.90, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[5].value, 114.14, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[6].value, -24.85, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[7].value,  31.00, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[8].value,  -0.65, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[9].value,  16.59, "%.2f" );
#if !RFC_MINIMAL
            /* Check matrix consistency */
            ASSERT( RFC_rfm_check( &ctx ) );
#endif /*!RFC_MINIMAL*/

        } while(0);
    }

    ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    while( mem_chain )
    {
        mem_chunk* next = mem_chain->next;
        free( mem_chain );
        mem_chain = next;
    }

    PASS();
}

#if !RFC_MINIMAL
TEST RFC_res_DIN45667( void )
{
/*
                              +
    8.5 ___________________________________________________
                                      +
    7.5 ___________________________________________________
                      +                   +
    6.5 ___________________________________________________
              +
    5.5 ___________________________________________________
          +                       +
    4.5 ___________________________________________________
                  +
    3.5 ___________________________________________________
                          +
    2.5 ___________________________________________________

    Slopes: 1, -2, 3, -4, 6, -4, 3, -1
    Sorted:  6,  3,  3,  1 ...  (raising)
            -4, -4, -2, -1      (falling)
    ==================================
    Counts: -4,  3,  2,  1
*/

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, -0.5 /* class_offset */,
                        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) )
    {
        RFC_VALUE_TYPE data[] = {4.9f, 6, 4, 7, 3, 9, 5, 8, 6.9f};
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );
    }
    else FAIL();

    ASSERT( ctx.residue_cnt == 8 );
    ASSERT( RFC_finalize( &ctx, RFC_RES_RP_DIN45667 ) );
    ASSERT( ctx.state == RFC_STATE_FINISHED );
    ASSERT( ctx.rp[0] == 0 );
    ASSERT( ctx.rp[1] == ctx.full_inc );
    ASSERT( ctx.rp[2] == ctx.full_inc );
    ASSERT( ctx.rp[3] == ctx.full_inc );
    ASSERT( ctx.rp[4] == ctx.full_inc );
    ASSERT( ctx.rp[5] == 0 );
    ASSERT( ctx.rp[6] == 0 );
    ASSERT( ctx.rp[7] == 0 );
    ASSERT( ctx.rp[8] == 0 );
    ASSERT( ctx.rp[9] == 0 );
    ASSERT( RFC_deinit( &ctx ) );

    PASS();
}

TEST RFC_res_repeated( void )
{
/*
                                                           
    8.5 ___________________________________________________
                      +               o                    
    7.5 ___________________________________________________
                                                           
    6.5 ___________________________________________________
              +               o                            
    5.5 ___________________________________________________
                                                           
    4.5 ___________________________________________________
                  +               o                        
    3.5 ___________________________________________________
          +               o                                
    2.5 ___________________________________________________

*/

#if RFC_TP_SUPPORT
    rfc_value_tuple_s tp[10];
#endif /*RFC_TP_SUPPORT*/
    double damage;
    double damage_5_3;
    double damage_7_2;

    damage_5_3 = pow( ( (5.0-3.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx;
    damage_7_2 = pow( ( (7.0-2.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx;
    damage = damage_5_3 * 2 + damage_7_2;

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, -0.5 /* class_offset */,
                        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) )
    {
        /*                       v--------v      */
        /*                          v--v         */
        RFC_VALUE_TYPE data[] = {2, 5, 3, 7};
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /*is_static*/ true ) );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
        ASSERT( RFC_dh_init( &ctx, RFC_SD_FULL_P2, /*dh*/ NULL, /*dh_cap*/ 0, /*is_static*/ true ) );
#endif /*RFC_DH_SUPPORT*/        
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );
    }
    else FAIL();

    ASSERT( ctx.state == RFC_STATE_BUSY_INTERIM );
    ASSERT( ctx.residue_cnt == 3 );
    ASSERT( ctx.damage == 0.0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 3 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT( ctx.tp[0].damage == 0.0 );
    ASSERT( ctx.tp[1].damage == 0.0 );
    ASSERT( ctx.tp[2].damage == 0.0 );
#endif /*RFC_DH_SUPPORT*/        
    ASSERT( ctx.residue[0].tp_pos == 1 );
    ASSERT( ctx.residue[1].tp_pos == 2 );
    ASSERT( ctx.residue[2].tp_pos == 3 );
#endif /*RFC_TP_SUPPORT*/

    ASSERT( RFC_finalize( &ctx, RFC_RES_REPEATED ) );
    ASSERT( ctx.state == RFC_STATE_FINISHED );
    ASSERT( ctx.residue_cnt == 0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 4 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
    ASSERT( ctx.tp[3].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT(       ctx.tp[0].damage == 0.0 );
    ASSERT( fabs( ctx.tp[1].damage / (damage_5_3*2) - 1 ) < 1e-10 );
    ASSERT(       ctx.tp[2].damage == 0.0 );
    ASSERT( fabs( ctx.tp[3].damage /  damage_7_2    - 1 ) < 1e-10 );
#endif /*RFC_DH_SUPPORT*/        
#endif /*RFC_TP_SUPPORT*/
    ASSERT( fabs( ctx.damage / damage - 1 ) < 1e-10 );
    ASSERT( RFC_deinit( &ctx ) );

    PASS();
}

TEST RFC_res_fullcycles( void )
{
/*
                                                           
    8.5 ___________________________________________________
                      +                                    
    7.5 ___________________________________________________
                                                           
    6.5 ___________________________________________________
              +                                            
    5.5 ___________________________________________________
                                                           
    4.5 ___________________________________________________
                  +                                        
    3.5 ___________________________________________________
          +                                                
    2.5 ___________________________________________________

*/

#if RFC_TP_SUPPORT
    rfc_value_tuple_s tp[10];
#endif /*RFC_TP_SUPPORT*/
    double damage;
    double damage_5_3 = pow( ( (5.0-3.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx;
    double damage_7_2 = pow( ( (7.0-2.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx;

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, -0.5 /* class_offset */,
                        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) )
    {
        /*                       v--------v      */
        /*                          v--v         */
        RFC_VALUE_TYPE data[] = {2, 5, 3, 7};
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /*is_static*/ true ) );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
        ASSERT( RFC_dh_init( &ctx, RFC_SD_HALF_23, /*dh*/ NULL, /*dh_cap*/ 0, /*is_static*/ true ) );
#endif /*RFC_DH_SUPPORT*/        
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );
    }
    else FAIL();

    ASSERT( ctx.state == RFC_STATE_BUSY_INTERIM );
    ASSERT( ctx.residue_cnt == 3 );
    ASSERT( ctx.damage == 0.0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 3 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT( ctx.tp[0].damage == 0.0 );
    ASSERT( ctx.tp[1].damage == 0.0 );
    ASSERT( ctx.tp[2].damage == 0.0 );
#endif /*RFC_DH_SUPPORT*/        
    ASSERT( ctx.residue[0].tp_pos == 1 );
    ASSERT( ctx.residue[1].tp_pos == 2 );
    ASSERT( ctx.residue[2].tp_pos == 3 );
#endif /*RFC_TP_SUPPORT*/

    ASSERT( RFC_finalize( &ctx, RFC_RES_FULLCYCLES ) );
    ASSERT( ctx.state == RFC_STATE_FINISHED );
    ASSERT( ctx.residue_cnt == 0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 4 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
    ASSERT( ctx.tp[3].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT( fabs( ctx.tp[0].damage / ( damage_7_2/2 ) - 1 ) < 1e-10 );
    ASSERT( fabs( ctx.tp[1].damage / ( damage_5_3/2 ) - 1 ) < 1e-10 );
    ASSERT( fabs( ctx.tp[2].damage / ( damage_5_3/2 ) - 1 ) < 1e-10 );
    ASSERT( fabs( ctx.tp[3].damage / ( damage_7_2/2 ) - 1 ) < 1e-10 );
#endif /*RFC_DH_SUPPORT*/        
#endif /*RFC_TP_SUPPORT*/
    damage = damage_7_2 + damage_5_3;
    ASSERT( fabs( ctx.damage / damage - 1 ) < 1e-10 );
    ASSERT( RFC_deinit( &ctx ) );

    PASS();
}

TEST RFC_res_halfcycles( void )
{
    /*

    8.5 ___________________________________________________
    +                                    
    7.5 ___________________________________________________

    6.5 ___________________________________________________
    +                                            
    5.5 ___________________________________________________

    4.5 ___________________________________________________
    +                                        
    3.5 ___________________________________________________
    +                                                
    2.5 ___________________________________________________

    */

#if RFC_TP_SUPPORT
    rfc_value_tuple_s tp[10];
#endif /*RFC_TP_SUPPORT*/
    double damage;
    double damage_5_3_half = pow( ( (5.0-3.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx / 2;
    double damage_7_2_half = pow( ( (7.0-2.0)/2 / ctx.wl_sx ), fabs(ctx.wl_k) ) / ctx.wl_nx / 2;

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, -0.5 /* class_offset */,
        1 /* hysteresis */, RFC_FLAGS_DEFAULT ) )
    {
        /*                       v--------v      */
        /*                          v--v         */
        RFC_VALUE_TYPE data[] = {2, 5, 3, 7};
#if RFC_TP_SUPPORT
        ASSERT( RFC_tp_init( &ctx, tp, NUMEL(tp), /*is_static*/ true ) );
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
        ASSERT( RFC_dh_init( &ctx, RFC_SD_HALF_23, /*dh*/ NULL, /*dh_cap*/ 0, /*is_static*/ true ) );
#endif /*RFC_DH_SUPPORT*/        
        ASSERT( RFC_feed( &ctx, data, NUMEL(data) ) );
    }
    else FAIL();

    ASSERT( ctx.state == RFC_STATE_BUSY_INTERIM );
    ASSERT( ctx.residue_cnt == 3 );
    ASSERT( ctx.damage == 0.0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 3 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT( ctx.tp[0].damage == 0.0 );
    ASSERT( ctx.tp[1].damage == 0.0 );
    ASSERT( ctx.tp[2].damage == 0.0 );
#endif /*RFC_DH_SUPPORT*/        
    ASSERT( ctx.residue[0].tp_pos == 1 );
    ASSERT( ctx.residue[1].tp_pos == 2 );
    ASSERT( ctx.residue[2].tp_pos == 3 );
#endif /*RFC_TP_SUPPORT*/

    ASSERT( RFC_finalize( &ctx, RFC_RES_HALFCYCLES ) );
    ASSERT( ctx.state == RFC_STATE_FINISHED );
    ASSERT( ctx.residue_cnt == 0 );
#if RFC_TP_SUPPORT
    ASSERT( ctx.tp_cnt == 4 );
    ASSERT( ctx.tp[0].tp_pos == 0 );
    ASSERT( ctx.tp[1].tp_pos == 0 );
    ASSERT( ctx.tp[2].tp_pos == 0 );
    ASSERT( ctx.tp[3].tp_pos == 0 );
#if RFC_DH_SUPPORT
    ASSERT( fabs( ctx.tp[0].damage / ( damage_7_2_half/2 ) - 1 ) < 1e-10 );  /* From residue */
    ASSERT( fabs( ctx.tp[1].damage / ( damage_5_3_half   ) - 1 ) < 1e-10 );  /* From regular counting */
    ASSERT( fabs( ctx.tp[2].damage / ( damage_5_3_half   ) - 1 ) < 1e-10 );  /* From regular counting */
    ASSERT( fabs( ctx.tp[3].damage / ( damage_7_2_half/2 ) - 1 ) < 1e-10 );  /* From residue */
#endif /*RFC_DH_SUPPORT*/        
#endif /*RFC_TP_SUPPORT*/
    damage = damage_7_2_half + damage_5_3_half * 2;
    ASSERT( fabs( ctx.damage / damage - 1 ) < 1e-10 );
    ASSERT( RFC_deinit( &ctx ) );

    PASS();
}
#endif /*!RFC_MINIMAL*/

#if RFC_AT_SUPPORT
double at_transform( rfc_ctx_s *rfc_ctx, double Sa, double Sm )
{
    double Sa_t;

    if( !RFC_at_transform( rfc_ctx, Sa, Sm, &Sa_t ) )
    {
        Sa_t = -1;
    }

    return Sa_t;
}

TEST RFC_at_test( void )
{
    bool ok = false;

    if( RFC_init( &ctx, 10 /* class_count */, 1 /* class_width */, 0 /* class_offset */,
                        1 /* hysteresis */, 0 /*flags*/ ) )
    {
        int i;
        double tol = 1e-10;
        double Sa[5], Sm[5];

        ok = true;

        ASSERT(
        RFC_at_init( &ctx, NULL /* Sa */, NULL /* Sm */, 0 /* count */, 0.3 /* M */, 
                           0.0 /* Sm_rig */, -1.0 /* R_rig */, true /* R_pinned */, false /* symmetric */ )
        );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/,  2.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/,  0.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/, -2.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.1 /*Sa*/,  9.0 /*Sm*/ ), 0.153636, 1e-5  /*tol*/ );  /* R =  89.0  */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/,  4.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/,  3.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  4.0 /*Sm*/ ), 2.836363, 1e-5  /*tol*/ );  /* R =  0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  3.0 /*Sm*/ ), 2.718181, 1e-5  /*tol*/ );  /* R =  0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  2.0 /*Sm*/ ), 2.6,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  3.0 /*Sm*/ ), 3.9,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  2.0 /*Sm*/ ), 3.6,      tol   /*tol*/ );  /* R = -0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  1.0 /*Sm*/ ), 2.3,      tol   /*tol*/ );  /* R = -0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  1.0 /*Sm*/ ), 3.3,      tol   /*tol*/ );  /* R = -0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  4.0 /*Sa*/,  1.0 /*Sm*/ ), 4.3,      tol   /*tol*/ );  /* R = -0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  0.0 /*Sm*/ ), 3.0,      tol   /*tol*/ );  /* R = -1.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -2.0 /*Sm*/ ), 1.4,      tol   /*tol*/ );  /* R = -Inf   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -9.0 /*Sm*/ ), 1.4,      tol   /*tol*/ );  /* R = -Inf   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.2 /*Sa*/, -9.0 /*Sm*/ ), 0.14,     tol   /*tol*/ );  /* R = -Inf   */

        ctx.at.R_rig = 0.6;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  1.0 /*Sm*/ ), 2.147928, 1e-5  /*tol*/ );  /* R = -0.5   */
        ctx.at.R_rig = 0.1;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  1.0 /*Sm*/ ), 2.488194, 1e-5  /*tol*/ );  /* R = -0.5   */
        ctx.at.R_rig = -4.0;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  1.0 /*Sm*/ ), 4.024390, 1e-5  /*tol*/ );  /* R = -0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  100.0 /*Sm*/ ), 5.6208425, 1e-5  /*tol*/ );  /* R ~ 0.942   */

        ASSERT(
        RFC_at_init( &ctx, NULL /* Sa */, NULL /* Sm */, 0 /* count */, 0.3 /* M */, 
                           0.0 /* Sm_rig */, -1.0 /* R_rig */, true /* R_pinned */, true /* symmetric */ )
        );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/,  2.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/,  0.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.0 /*Sa*/, -2.0 /*Sm*/ ), 0.0,      0.0   /*tol*/ );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  0.1 /*Sa*/,  9.0 /*Sm*/ ), 0.153636, 1e-5  /*tol*/ );  /* R =  89.0  */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/,  4.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/,  3.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  4.0 /*Sm*/ ), 2.836363, 1e-5  /*tol*/ );  /* R =  0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  3.0 /*Sm*/ ), 2.718181, 1e-5  /*tol*/ );  /* R =  0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  2.0 /*Sm*/ ), 2.6,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  3.0 /*Sm*/ ), 3.9,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  2.0 /*Sm*/ ), 3.6,      tol   /*tol*/ );  /* R = -0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/,  1.0 /*Sm*/ ), 2.3,      tol   /*tol*/ );  /* R = -0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  1.0 /*Sm*/ ), 3.3,      tol   /*tol*/ );  /* R = -0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  4.0 /*Sa*/,  1.0 /*Sm*/ ), 4.3,      tol   /*tol*/ );  /* R = -0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/,  0.0 /*Sm*/ ), 3.0,      tol   /*tol*/ );  /* R = -1.0   */

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/, -4.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/, -3.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -4.0 /*Sm*/ ), 2.836363, 1e-5  /*tol*/ );  /* R =  0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -3.0 /*Sm*/ ), 2.718181, 1e-5  /*tol*/ );  /* R =  0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -2.0 /*Sm*/ ), 2.6,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -3.0 /*Sm*/ ), 3.9,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -2.0 /*Sm*/ ), 3.6,      tol   /*tol*/ );  /* R = -0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -1.0 /*Sm*/ ), 2.3,      tol   /*tol*/ );  /* R = -0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -1.0 /*Sm*/ ), 3.3,      tol   /*tol*/ );  /* R = -0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  4.0 /*Sa*/, -1.0 /*Sm*/ ), 4.3,      tol   /*tol*/ );  /* R = -0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -0.0 /*Sm*/ ), 3.0,      tol   /*tol*/ );  /* R = -1.0   */

        for( i = 0; i < 5; i++ )
        {
            Sa[i] = ctx.at.Sa[i] * 333;
            Sm[i] = ctx.at.Sm[i] * 333;
        }

        ASSERT(
        RFC_at_init( &ctx, Sa /* Sa */, Sm /* Sm */, 5 /* count */, 0.3 /* M */, 
                           0.0 /* Sm_rig */, -1.0 /* R_rig */, true /* R_pinned */, false /* symmetric */ )
        );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/, -4.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  1.0 /*Sa*/, -3.0 /*Sm*/ ), 1.536363, 1e-5  /*tol*/ );  /* R =  0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -4.0 /*Sm*/ ), 2.836363, 1e-5  /*tol*/ );  /* R =  0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -3.0 /*Sm*/ ), 2.718181, 1e-5  /*tol*/ );  /* R =  0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -2.0 /*Sm*/ ), 2.6,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -3.0 /*Sm*/ ), 3.9,      tol   /*tol*/ );  /* R =  0.0   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -2.0 /*Sm*/ ), 3.6,      tol   /*tol*/ );  /* R = -0.2   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  2.0 /*Sa*/, -1.0 /*Sm*/ ), 2.3,      tol   /*tol*/ );  /* R = -0.333 */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -1.0 /*Sm*/ ), 3.3,      tol   /*tol*/ );  /* R = -0.5   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  4.0 /*Sa*/, -1.0 /*Sm*/ ), 4.3,      tol   /*tol*/ );  /* R = -0.6   */
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  3.0 /*Sa*/, -0.0 /*Sm*/ ), 3.0,      tol   /*tol*/ );  /* R = -1.0   */


        ASSERT(
        RFC_at_init( &ctx, NULL /* Sa */, NULL /* Sm */, 0 /* count */, 0.3 /* M */, 
                           50.0 /* Sm_rig */, 0.0 /* R_rig */, false /* R_pinned */, false /* symmetric */ )
        );

        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  100.0 /*Sa*/, 0.0 /*Sm*/ ), 85.0, tol  /*tol*/ );
        ctx.at.Sm_rig = 200.0;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  100.0 /*Sa*/, 0.0 /*Sm*/ ), 65.088757, 1e-5  /*tol*/ );
        ctx.at.Sm_rig = 400.0;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx,  100.0 /*Sa*/, 50.0 /*Sm*/ ), 74.85207, 1e-5  /*tol*/ );

        for( i = 0; i < 3; i++ )
        {
            Sa[i] = ctx.at.Sa[i] * 333;
            Sm[i] = ctx.at.Sm[i] * 333;
        }
        ctx.at.Sm_rig = 400.0;
        GREATEST_ASSERT_IN_RANGE( at_transform( &ctx, 100.0 /*Sa*/, 50.0 /*Sm*/ ), 74.85207, 1e-5  /*tol*/ );

    }

    ASSERT( RFC_deinit( &ctx ) );

    if( ok )
    {
        PASS();
    }
    else
    {
        FAIL();
    }
}
#endif /*RFC_AT_SUPPORT*/

TEST RFC_CPP_wrapper_easy( void )
{
    extern bool wrapper_test_easy( void );
    ASSERT( wrapper_test_easy() );
    PASS();
}

TEST RFC_CPP_wrapper_advanced( void )
{
    extern bool wrapper_test_advanced( void );
    ASSERT( wrapper_test_advanced() );
    PASS();
}

#if !RFC_MINIMAL
TEST RFC_wl_math( void )
{
    double sd =  1e3;
    double nd =  1e7;
    double k  = -5;
    double k2 = -9;
    double sx =  300;
    double nx =  pow( sx/sd, k2 ) * nd;
    double s0 =  500;
    double n0 =  pow( s0/sx, k ) * nx;
    double x;

    x = ( log(nx) - log(nd) ) / ( log(sx) - log(sd) );
    ASSERT_IN_RANGE( k2, x, 1e-3 );

    x = ( log(n0) - log(nx) ) / ( log(s0) - log(sx) );
    ASSERT_IN_RANGE( k, x, 1e-3 );

    ASSERT( RFC_wl_calc_k2( &ctx, s0, n0, k, sx, nx, &x, sd, nd ) );
    ASSERT_IN_RANGE( k2, x, 1e-3 );

    ASSERT( RFC_wl_calc_sx( &ctx, s0, n0, k, &x, nx, k2, sd, nd ) );
    ASSERT_IN_RANGE( sx, x, 1e-3 );

    ASSERT( RFC_wl_calc_sd( &ctx, s0, n0, k, sx, nx, k2, &x, nd ) );
    ASSERT_IN_RANGE( sd, x, 1e-3 );

    ASSERT( RFC_wl_calc_sa( &ctx, sx, nx, k2, nd, &x ) );
    ASSERT_IN_RANGE( sd, x, 1e-3 );

    ASSERT( RFC_wl_calc_n( &ctx, sx, nx, k2, sd, &x ) );
    ASSERT_IN_RANGE( nd, x, 1e-3 );

    PASS();
}
#endif /*!RFC_MINIMAL*/


#if !RFC_MINIMAL
TEST RFC_miner_consequent( void )
{
    /* Data from [6] table 3.2-6 */
    /* =================================================================================================================== */
    /* Amplitude-counts histogram */
    double          Sa_rel[]        = { 0.000,  0.125,   0.275,   0.425,  0.575,  0.725,  0.850,  0.950,  1.000 };
    RFC_counts_type Sa_counts[]     = { 0,      605000,  280000,  92000,  20000,  2720,   280,    16,     2 };
    /* Various representations for the histogram defined above */
    double          Sa_hat[]        = { 100, 105, 110, 115, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800 };
    double          A_expected[]    = { 89199.590, 24445.830, 14414.850, 6980.954, 2089.658, 556.181, 253.551, 219.482, 
                                        152.775, 144.658, 133.296, 129.536, 129.245, 128.810, 128.205, 127.398 };
    /* Woehler curve */
    const double SD = 100.0;
    const double ND = 1e6;
    const double k  = 4;
    /* =================================================================================================================== */

    unsigned        class_count     = NUMEL(Sa_counts);
    RFC_value_type  class_offset    = 0.0,
                    hysteresis      = 0.0;
    double          D_mk;
    int             i;

    /* Check each representation */
    for( i = 0; i < NUMEL(Sa_hat); i++ )
    {
        RFC_value_type  class_width = Sa_hat[i] * 2.0 / ( class_count - 1 );
        double          Sa[10];
        int             j;
        double          N_bar;
        double          A;
        double          h_sum;

        N_bar = pow( SD / Sa_hat[i], k ) * ND;

        ASSERT( NUMEL(Sa) >= class_count );

        h_sum = 0;
        for( j = 0; j < (int)class_count; j++ )
        {
            Sa[j]  = Sa_rel[j] * Sa_hat[i];
            h_sum += Sa_counts[j];
        }

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
        ASSERT( RFC_wl_init_original( &ctx, SD, ND, k ) );

        ctx.full_inc = 1;
        ASSERT( RFC_damage_from_rp( &ctx, Sa_counts, Sa, &D_mk, RFC_RP_DAMAGE_CALC_TYPE_CONSEQUENT /*rp_calc_type*/ ) );

        /* A is the difference from variable-amplitude to constant-amplitude fatigue life for Sa_hat.
           (Sa_hat is the largest given Sa in histogram) */
        A = h_sum / D_mk / N_bar;

        /* Check with the proposed values given in [6] */
        GREATEST_ASSERT_IN_RANGE( A_expected[i], A, 0.1 );

        ASSERT( RFC_deinit( &ctx ) );
    }
    PASS();
}

TEST RFC_miner_consequent2( void )
{
    RFC_VALUE_TYPE      data[10000];
    size_t              data_len            =  NUMEL( data );
    RFC_VALUE_TYPE      x_max;
    RFC_VALUE_TYPE      x_min;
    unsigned            class_count         =  100;
    RFC_VALUE_TYPE      class_width;
    RFC_VALUE_TYPE      class_offset;
    RFC_VALUE_TYPE      hysteresis;
    double              D_elem, D_orig, D_mod, D_con;
    double              rp_hist[100];
    double              k, k2;
    size_t              i, j;
    size_t              repeats;

    ASSERT( NUMEL(rp_hist) == class_count );

#include "long_series.c"

    for( i = 0; i < data_len; i++ )
    {
        double value = data_export[i];
        data[i] = value;
        if( !i )
        {
            x_max = x_min = value;
        }
        else
        {
            if( value > x_max ) x_max = value;
            if( value < x_min ) x_min = value;
        }
    }

    class_width  =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    class_offset =  x_min - class_width / 2;
    hysteresis   =  class_width;

    ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
    ASSERT( RFC_wl_init_original( &ctx, /*sd */ x_max * 0.075, /* nd */ 1e5, /* k */ -1.5 ) );

    k  = -fabs( ctx.wl_k );
    k2 = -( 2 * fabs( ctx.wl_k ) - 1 );

    for( i = 0; ctx.internal.wl.D < 1.0; i++ )
    {
        ASSERT( RFC_feed_scaled( &ctx, data, /* count */ data_len, /* factor */ 1.0 ) );

        for( j = 0; j < 0; j++ )
        {
            ASSERT( RFC_feed_scaled( &ctx, data, /* count */ data_len, /* factor */ 0.02 ) );
        }
    }

    repeats = i;

    ASSERT( RFC_damage_from_rp( &ctx, /* counts */ NULL, /* sa */ NULL, &D_orig, RFC_RP_DAMAGE_CALC_TYPE_DEFAULT ) );
    ASSERT( RFC_damage_from_rp( &ctx, /* counts */ NULL, /* sa */ NULL, &D_elem, RFC_RP_DAMAGE_CALC_TYPE_ELEMENTAR ) );
    ctx.wl_k2 = k2;
    ASSERT( RFC_damage_from_rp( &ctx, /* counts */ NULL, /* sa */ NULL, &D_mod,  RFC_RP_DAMAGE_CALC_TYPE_MODIFIED ) );
    ctx.wl_k2 = k;
    ASSERT( RFC_damage_from_rp( &ctx, /* counts */ NULL, /* sa */ NULL, &D_con,  RFC_RP_DAMAGE_CALC_TYPE_CONSEQUENT ) );

    ASSERT( fabs( D_con / ctx.internal.wl.D - 1 ) < 1e-3 );

    for( j = 0; j < (int)class_count; j++ )
    {
        rp_hist[j] = (double)ctx.rp[j] / ctx.full_inc;
    }

    ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

    if(0)
    {
        FILE*   file = NULL;
        int     i;
        double  N = 0.0;

        setlocale( LC_ALL, "" );

        file = fopen( "miner_consequent_results.txt", "wt" );
        ASSERT( file );
        fprintf( file, "Class parameters:\n" );
        fprintf( file, "Class count: \t%d\n", (int)ctx.class_count );
        fprintf( file, "Class width:  \t%.5f\n", ctx.class_width );
        fprintf( file, "Class offset:  \t%.5f\n", ctx.class_offset );
        fprintf( file, "\nWoehler parameters:\n" );
        fprintf( file, "Sx: \t%g\n", ctx.wl_sx );
        fprintf( file, "Nx: \t%g\n", ctx.wl_nx );
        fprintf( file, "k: \t%g\n", k );
        fprintf( file, "Sd: \t%g\n", ctx.wl_sd );
        fprintf( file, "Nd: \t%g\n", ctx.wl_nd );
        fprintf( file, "k2: \t%g\n", k2 );
        fprintf( file, "Omission: \t%g\n", ctx.wl_omission );
        fprintf( file, "\nDamage:\n" );
        fprintf( file, "Miner original: \t%g\n", D_orig );
        fprintf( file, "Miner elementary: \t%g\n", D_elem );
        fprintf( file, "Miner modified: \t%g\n", D_mod );
        fprintf( file, "Miner consequent: \t%g\n", D_con );
        fprintf( file, "Miner live: \t%g\n", ctx.internal.wl.D );
        fprintf( file, "\nRepeats: \t%lu\n", (unsigned long)repeats );
        fprintf( file, "\nRP Histogram:\n" );
        fprintf( file, "#\tSa\tn\tN\tWL_orig\tWL_elem\tWL_mod\n" );

        for( i = (int)ctx.class_count - 1; i >= 0; i-- )
        {
            double Sa = (double)ctx.class_width * i / 2;
            double n = rp_hist[i];
            double x;

            N += n;

            if( N == 0.0 ) continue;

            fprintf( file, "%d\t", i );
            fprintf( file, "%g\t", Sa );
            fprintf( file, "%g\t", n );
            fprintf( file, "%g\t", N );

            if( Sa < ctx.wl_sd )
            {
                x = 1e30; /*DBL_MAX*/
            }
            else
            {
                ASSERT( RFC_wl_calc_n( &ctx, ctx.wl_sx, ctx.wl_nx, k, Sa, &x ) );
            }
            fprintf( file, "%g\t", x );

            if( Sa > 0.0 )
            {
                ASSERT( RFC_wl_calc_n( &ctx, ctx.wl_sx, ctx.wl_nx, k, Sa, &x ) );
            }
            else
            {
                x = 1e30; /*DBL_MAX*/
            }
            fprintf( file, "%g\t", x );

            if( Sa > ctx.wl_sx )
            {
                ASSERT( RFC_wl_calc_n( &ctx, ctx.wl_sx, ctx.wl_nx, k, Sa, &x ) );
            }
            else if( Sa > 0.0 )
            {
                ASSERT( RFC_wl_calc_n( &ctx, ctx.wl_sx, ctx.wl_nx, k2, Sa, &x ) );
            }
            else
            {
                x = 1e30; /*DBL_MAX*/
            }
            fprintf( file, "%g\t", x );
            fprintf( file, "\n" );
        }
        fclose( file );
    }

    ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}
#endif /*!RFC_MINIMAL*/


/* local suite (greatest) */
SUITE( RFC_TEST_SUITE )
{
    /* Test rainflow counting */
    RUN_TEST1( RFC_empty, 1 );
    RUN_TEST1( RFC_cycle_up, 1 );
    RUN_TEST1( RFC_cycle_down, 1 );
    RUN_TEST1( RFC_small_example, 1 );
    RUN_TEST1( RFC_long_series, 1 );
    RUN_TEST1( RFC_empty, 0 );
    RUN_TEST1( RFC_cycle_up, 0 );
    RUN_TEST1( RFC_cycle_down, 0 );
    RUN_TEST1( RFC_small_example, 0 );
    RUN_TEST1( RFC_long_series, 0 );
    /* Test C++ Wrapper */
    RUN_TEST( RFC_CPP_wrapper_easy );
    RUN_TEST( RFC_CPP_wrapper_advanced );
#if !RFC_MINIMAL
    /* Residual methods */
    RUN_TEST( RFC_res_DIN45667 );
    RUN_TEST( RFC_res_repeated );
    RUN_TEST( RFC_res_fullcycles );
    RUN_TEST( RFC_res_halfcycles );
    RUN_TEST( RFC_wl_math );
    /* "Miner consequent" approach */
    RUN_TEST( RFC_miner_consequent );
    RUN_TEST( RFC_miner_consequent2 );
#endif /*!RFC_MINIMAL*/
#if RFC_TP_SUPPORT
    /* Test turning points */
    RUN_TEST1( RFC_test_turning_points, 1 );
    RUN_TEST1( RFC_tp_prune_test, 1 );
    RUN_TEST1( RFC_test_turning_points, 0 );
    RUN_TEST1( RFC_tp_prune_test, 0 );
#endif /*RFC_TP_SUPPORT*/
#if RFC_AT_SUPPORT
    /* Test amplitude transformation */
    RUN_TEST( RFC_at_test );
#endif /*RFC_AT_SUPPORT*/
}


/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();


int main( int argc, char *argv[] )
{
    if( argc > 1 )
    {
        FILE* file;

        file = fopen(  argv[1], "rt" );
        if( file )
        {
            fclose( file );
            long_series_file = argv[1];
        }
    }

    if( !long_series_file )
    {
        long_series_file = "long_series.txt";
    }

    GREATEST_MAIN_BEGIN();      /* init & parse command-line args */
    RUN_SUITE( RFC_TEST_SUITE );
    GREATEST_MAIN_END();        /* display results */

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    while( mem_chain )
    {
        mem_chunk* next = mem_chain->next;
        free( mem_chain );
        mem_chain = next;
    }
}
