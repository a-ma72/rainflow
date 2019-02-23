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
 *    Test suite
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

#define GREATEST_FPRINTF fprintf
#define RFC_VALUE_TYPE   double

#include "../rainflow.h"
#include "../greatest/greatest.h"
#include <locale.h>
#include <math.h>
#include <float.h>

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define NUMEL(x) (sizeof(x)/sizeof((x)[0]))

static
struct buffer
{
    void *ptr;
    struct buffer *next;
} buffers = {0};


typedef struct mem_chunk
{
    size_t             size, 
                       count;
    struct  mem_chunk *next;
    RFC_VALUE_TYPE     data[1];
} mem_chunk;

      rfc_ctx_s   ctx              = { sizeof(ctx) };   /* module shared rainflow context */
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

static 
struct buffer* add_buffer( void* ptr )
{
    struct buffer *new_buffer = (struct buffer*)calloc( 1, sizeof( struct buffer ) );
    struct buffer *buffer_ptr = &buffers;

    while( buffer_ptr->next )
    {
        buffer_ptr = buffer_ptr->next;
    }

    new_buffer->ptr  = ptr;
    buffer_ptr->next = new_buffer;

    return buffer_ptr;
}

static
void strip_buffer( struct buffer* buffer )
{
    struct buffer* buffer_next;

    if( !buffer ) buffer = &buffers;

    buffer_next  = buffer->next;
    buffer->next = NULL;

    while( buffer_next )
    {
        buffer = buffer_next;
        buffer_next = buffer->next;
        if( buffer->ptr )
        {
            free( buffer->ptr );
        }
        free( buffer );
    }
}




rfc_counts_t rfm_peek( rfc_ctx_s *rfc_ctx, int from, int to )
{
    from = (int)( ( (double)from - rfc_ctx->class_offset ) / rfc_ctx->class_width );
    to   = (int)( ( (double)to   - rfc_ctx->class_offset ) / rfc_ctx->class_width );
    
    return rfc_ctx->rfm[ from * rfc_ctx->class_count + to ];
}




TEST RFC_empty( int ccnt )
{
    RFC_VALUE_TYPE      x_max           =  1;
    RFC_VALUE_TYPE      x_min           = -1;
    unsigned            class_count     =  ccnt ? 100 : 0;
    RFC_VALUE_TYPE      class_width     =  (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
    RFC_VALUE_TYPE      class_offset    =  x_min - class_width / 2;
    RFC_VALUE_TYPE      hysteresis      =  class_width;
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {0};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
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
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {1,3,2,4};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 3, 2 ), 1 * ctx.full_inc );
            ASSERT_EQ( ctx.residue_cnt, 2 );
            ASSERT_EQ( ctx.residue[0].value, 1.0 );
            ASSERT_EQ( ctx.residue[1].value, 4.0 );
            ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
            ASSERT_EQ( ctx.residue[0].pos, 1 );
            ASSERT_EQ( ctx.residue[1].pos, 4 );
        }
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
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {4,2,3,1};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 1.0 );
            ASSERT_EQ( rfm_peek( &ctx, 2, 3 ), 1 * ctx.full_inc );
            ASSERT_EQ( ctx.residue_cnt, 2 );
            ASSERT_EQ( ctx.residue[0].value, 4.0 );
            ASSERT_EQ( ctx.residue[1].value, 1.0 );
            ASSERT_EQ( ctx.state, RFC_STATE_FINISHED );
            ASSERT_EQ( ctx.residue[0].pos, 1 );
            ASSERT_EQ( ctx.residue[1].pos, 4 );
        }
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
    size_t              i;

    do
    {
        RFC_VALUE_TYPE data[] = {2,5,3,6,2,4,1,6,1,4,1,5,3,6,3,6,1,5,2};
        RFC_VALUE_TYPE sum = 0.0;

        ASSERT( RFC_init( &ctx, class_count, class_width, class_offset, hysteresis, RFC_FLAGS_DEFAULT ) );
        ASSERT( RFC_feed( &ctx, data, /* count */ NUMEL( data ) ) );
        ASSERT( RFC_finalize( &ctx, /* residual_method */ RFC_RES_NONE ) );

        for( i = 0; i < class_count * class_count; i++ )
        {
            sum += ctx.rfm[i] / ctx.full_inc;
        }

        if( class_count )
        {
            ASSERT_EQ( sum, 7.0 );
            ASSERT_EQ( rfm_peek( &ctx, 5, 3 ), 2 * ctx.full_inc );
            ASSERT_EQ( rfm_peek( &ctx, 6, 3 ), 1 * ctx.full_inc );
            ASSERT_EQ( rfm_peek( &ctx, 1, 4 ), 1 * ctx.full_inc );
            ASSERT_EQ( rfm_peek( &ctx, 2, 4 ), 1 * ctx.full_inc );
            ASSERT_EQ( rfm_peek( &ctx, 1, 6 ), 2 * ctx.full_inc );
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
        GREATEST_FPRINTF( GREATEST_STDOUT, "Class count (%d): ", class_count );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) && ( value > 0.0 ) )
            {
                class_count = (unsigned)( value + 0.5 );
            }
        }

        hysteresis = class_width = (RFC_VALUE_TYPE)ROUND( 100 * (x_max - x_min) / (class_count - 1) ) / 100;
        GREATEST_FPRINTF( GREATEST_STDOUT, "Class width (%g): ", class_width );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) && ( value > 0.0 ) )
            {
                hysteresis = class_width = value;
            }
        }

        class_offset  =  x_min - class_width / 2;
        GREATEST_FPRINTF( GREATEST_STDOUT, "Class offset (%g): ", class_offset );
        if( fgets( buf, sizeof(buf), stdin ) != NULL )
        {
            if( ( 1 == sscanf( buf, "%lf %n", &value, &len ) ) && ( strlen(buf) == len ) )
            {
                class_offset = value;
            }
        }
        GREATEST_FPRINTF( GREATEST_STDOUT, "\n" );
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
            /* Check residue */
            ASSERT_EQ( ctx.residue_cnt, 10 );
            ASSERT_EQ_FMT( ctx.residue[0].value,   0.54, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[1].value,   2.37, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[2].value,  -0.45, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[3].value,  17.04, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[4].value, -50.90, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[5].value, 114.14, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[6].value, -24.85, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[7].value,  31.00, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[8].value,  -0.65, "%.2f" );
            ASSERT_EQ_FMT( ctx.residue[9].value,  16.59, "%.2f" );
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

    /* Test C++ Wrapper */
    GREATEST_SUITE_EXTERN( RFC_WRAPPER_SUITE_SIMPLE );
    RUN_SUITE( RFC_WRAPPER_SUITE_SIMPLE );
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

    strip_buffer( NULL );
}
