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

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define NUMEL(x) (sizeof(x)/sizeof((x)[0]))

rfc_ctx_s ctx = { sizeof(ctx) };

double rfm_peek( rfc_ctx_s *rfc_ctx, int from, int to )
{
    return (double)rfc_ctx->matrix[ (from-1)*rfc_ctx->class_count + (to-1)] / rfc_ctx->full_inc;
}


TEST RFC_empty(void)
{
    RFC_VALUE_TYPE      x_max           =  1;
    RFC_VALUE_TYPE      x_min           = -1;
    unsigned            class_count     =  100;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width;
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {0};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ 0 ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.matrix[i];
        }

        ASSERT_EQ( sum, 0.0 );
        ASSERT_EQ( ctx.residue_cnt, 0 );
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}


TEST RFC_cycle_up(void)
{
    RFC_VALUE_TYPE      x_max           =  4;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  4;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width * 0.99;
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {1,3,2,4};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.matrix[i] / ctx.full_inc;
        }

        ASSERT_EQ( sum, 1.0 );
        ASSERT_EQ( rfm_peek( &ctx, 3, 2 ), 1.0 );
        ASSERT_EQ( ctx.residue_cnt, 2 );
        ASSERT_EQ( ctx.residue[0].value, 1.0 );
        ASSERT_EQ( ctx.residue[1].value, 4.0 );
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}


TEST RFC_cycle_down(void)
{
    RFC_VALUE_TYPE      x_max           =  4;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  4;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width * 0.99;
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {4,2,3,1};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.matrix[i] / ctx.full_inc;
        }

        ASSERT_EQ( sum, 1.0 );
        ASSERT_EQ( rfm_peek( &ctx, 2, 3 ), 1.0 );
        ASSERT_EQ( ctx.residue_cnt, 2 );
        ASSERT_EQ( ctx.residue[0].value, 4.0 );
        ASSERT_EQ( ctx.residue[1].value, 1.0 );
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}


TEST RFC_small_example(void)
{
    RFC_VALUE_TYPE      x_max           =  6;
    RFC_VALUE_TYPE      x_min           =  1;
    unsigned            class_count     =  (unsigned)x_max;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width;
    size_t          i;

    do
    {
        RFC_VALUE_TYPE data[] = {2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.matrix[i] / ctx.full_inc;
        }

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
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}


TEST RFC_long_series(void)
{
    FILE*               file            =  NULL;
    RFC_VALUE_TYPE      data[10000];
    size_t              data_len        =  NUMEL( data );
    RFC_VALUE_TYPE      x_max;
    RFC_VALUE_TYPE      x_min;
    unsigned            class_count     =  100;
    RFC_VALUE_TYPE      class_width;
    RFC_VALUE_TYPE      class_offset;
    RFC_VALUE_TYPE      hysteresis;
    size_t              i;

    file = fopen( "long_series.csv", "rt" );
    ASSERT( file );
    for( i = 0; i < data_len; i++ )
    {
        double value;
        fscanf( file, "%lf\n", &value );
        data[i] = value;
        if( !i )
        {
            x_max = x_min = value;
        }
        else
        {
            if( value > x_max ) x_max = data[i];
            if( value < x_min ) x_min = data[i];
        }
    }
    fclose( file );

    class_width   =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    class_offset  =  x_min - class_width / 2;
    hysteresis    =  class_width;

    do
    {
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ data_len ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.matrix[i] / ctx.full_inc;
        }

        ASSERT_EQ( sum, 602.0 );
        GREATEST_ASSERT_IN_RANGE( ctx.pseudo_damage, 4.8703e-16, 0.00005e-16 );
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
        ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
    } while(0);

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }

    PASS();
}



/* local suite (greatest) */
SUITE( RFC_TEST_SUITE )
{
    RUN_TEST( RFC_empty );
    RUN_TEST( RFC_cycle_up );
    RUN_TEST( RFC_cycle_down );
    RUN_TEST( RFC_small_example );
    RUN_TEST( RFC_long_series );
}


/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();


int main( int argc, char *argv[] )
{
    GREATEST_MAIN_BEGIN();      /* init & parse command-line args */
    RUN_SUITE( RFC_TEST_SUITE );
    GREATEST_MAIN_END();        /* display results */

    if( ctx.state != RFC_STATE_INIT0 )
    {
        RFC_deinit( &ctx );
    }
}
