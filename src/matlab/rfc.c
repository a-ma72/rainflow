#include "rainflow.h"

#include <ctype.h>
#include <string.h>
#include <mex.h>

#define MAT_OFFS( i, j )  ( (i) * class_count + (j) )
#define CALLOC calloc
#define REALLOC realloc
#define FREE free

static void * mem_alloc( void *ptr, size_t num, size_t size, int aim );

#if !RFC_MINIMAL
#define RFC_MEX_USAGE \
"\nUsage:\n"\
"[pd,re,rm,rp,lc,tp,dh] = rfc( 'rfc', data, class_count, class_width, class_offset, hysteresis, residual_method, enforce_margin, use_hcm )\n"\
"    pd = Pseudo damage\n"\
"    re = Residue\n"\
"    rm = Rainflow matrix (from/to)\n"\
"    rp = Range pair counts\n"\
"    lc = Level crossings\n"\
"    tp = Turning points\n"\
"    dh = Damage history\n"\
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



#if RFC_DEBUG_FLAGS
static
int rfc_vfprintf_fcn( void *ctx, FILE* stream, const char *fmt, va_list arg )
{
    int length;

    /* Get the buffer size needed */
    length = vsnprintf( NULL, 0, fmt, arg );  
    if( length > 0 )
    {
        char *buffer = CALLOC( ++length, 1 );

        if( buffer )
        {
            buffer[length-1] = 0;
            vsnprintf( buffer, length, fmt, arg );
            mexPrintf( "%s", buffer );

            FREE( buffer );
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
    if( nrhs != 10 + 1 )
    {
        mexErrMsgTxt( "Function needs exact 10 arguments!" );
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
        const mxArray  *mxUseASTM        = prhs[8];
        const mxArray  *mxSpreadDamage   = prhs[9];
        const mxArray  *mxAutoResize     = prhs[10];
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
        int             use_astm         = (int)mxGetScalar( mxUseASTM );
        int             spread_damage    = (int)mxGetScalar( mxSpreadDamage );
        int             auto_resize      = (int)mxGetScalar( mxAutoResize );
#else /*RFC_MINIMAL*/
        int             residual_method  = RFC_RES_NONE;
#endif /*!RFC_MINIMAL*/
        size_t          i;
        bool            ok;

        mxAssert( residual_method >= 0 && residual_method < RFC_RES_COUNT, 
                  "Invalid residual method!" );

#if !RFC_MINIMAL
#if RFC_HCM_SUPPORT
        mxAssert( use_hcm == 0 || use_hcm == 1,
                  "Invalid HCM flag, use 0 or 1!" );
#else /*!RFC_HCM_SUPPORT*/
        mxAssert( use_hcm == 0,
                  "HCM not supported!" );
#endif /*RFC_HCM_SUPPORT*/

#if RFC_ASTM_SUPPORT
        mxAssert( use_astm == 0 || use_astm == 1,
                  "Invalid ASTM flag, use 0 or 1!" );
#if RFC_HCM_SUPPORT
        mxAssert( use_astm == 0 || use_hcm == 0,
                  "Invalid ASTM flag, use either HCM or ASTM!" );
#endif /*RFC_HCM_SUPPORT*/
#else /*!RFC_ASTM_SUPPORT*/
        mxAssert( use_astm == 0,
                  "ASTM not supported!" );
#endif /*RFC_ASTM_SUPPORT*/


#if RFC_DH_SUPPORT
        mxAssert( spread_damage >= RFC_SD_NONE && spread_damage < RFC_SD_COUNT,
                  "Invalid spread damage method!" );
#else /*!RFC_DH_SUPPORT*/
        if( spread_damage != 0 )
        {
            mexErrMsgTxt( "Invalid spread damage method, only 0 accepted!" );
        }
#endif /*RFC_DH_SUPPORT*/

#if RFC_AR_SUPPORT
        mxAssert( auto_resize == 0 || auto_resize == 1,
                  "Invalid auto resize flag, use 0 or 1!" );
#else
        if( auto_resize != 0 )
        {
            mexErrMsgTxt( "Invalid auto resize flag, only 0 accepted!" );
        }
#endif /*RFC_AR_SUPPORT*/
#endif /*!RFC_MINIMAL*/

        ok = RFC_init( &rfc_ctx, 
                       class_count, (rfc_value_t)class_width, (rfc_value_t)class_offset, 
                       (rfc_value_t)hysteresis, RFC_FLAGS_DEFAULT );
        if( !ok )
        {
            RFC_deinit( &rfc_ctx );
            mexErrMsgTxt( "Error during initialization!" );
        }

#if RFC_DEBUG_FLAGS
        rfc_ctx.debug_vfprintf_fcn = rfc_vfprintf_fcn;
#endif /*RFC_DEBUG_FLAGS*/

#if RFC_TP_SUPPORT
        ok = RFC_tp_init( &rfc_ctx, /*tp*/ NULL, /*tp_cap*/ 128, /* is_static */ false );

        if( !ok )
        {
            RFC_deinit( &rfc_ctx );
            mexErrMsgTxt( "Error during initialization (tp)!" );
        }
#endif /*RFC_TP_SUPPORT*/

        /* Cast values from double type to rfc_value_t */ 
        if( sizeof( rfc_value_t ) != sizeof(double) && data_len )  /* maybe unsafe! */
        {
            buffer = (rfc_value_t *)mem_alloc( NULL, data_len, 
                                               sizeof(rfc_value_t), RFC_MEM_AIM_TEMP );

            if( !buffer )
            {
                RFC_deinit( &rfc_ctx );
                mexErrMsgTxt( "Error during initialization (memory)!" );
            }

            for( i = 0; i < data_len; i++ )
            {
                buffer[i] = (rfc_value_t)data[i];
            }
        }
        else buffer = (rfc_value_t*)data;

#if RFC_DH_SUPPORT
        if( spread_damage >= RFC_SD_TRANSIENT_23 )
        {
            if( !RFC_dh_init( &rfc_ctx, spread_damage, /*dh*/ NULL, /*dh_cap*/ 1, /*is_static*/ false ) )
            {
                ok = false;
            }
        }
        else
        {
            if( !RFC_dh_init( &rfc_ctx, spread_damage, /*dh*/ NULL, /*dh_cap*/ 0, /*is_static*/ true ) )
            {
                ok = false;
            }
        }
        
#endif /*RFC_DH_SUPPORT*/

#if RFC_AR_SUPPORT
        if( auto_resize )
        {
            RFC_flags_set( &rfc_ctx, RFC_FLAGS_AUTORESIZE, /* stack */ 0, /* overwrite */ false );
        }
#endif /*RFC_AR_SUPPORT*/

        /* Rainflow counting */

#if !RFC_MINIMAL
        /* Setup */
        rfc_ctx.internal.flags  |= enforce_margin ? RFC_FLAGS_ENFORCE_MARGIN : 0;
#endif /*!RFC_MINIMAL*/

#if RFC_HCM_SUPPORT && RFC_ASTM_SUPPORT
             if( use_hcm )  rfc_ctx.counting_method = RFC_COUNTING_METHOD_HCM;
        else if( use_astm ) rfc_ctx.counting_method = RFC_COUNTING_METHOD_ASTM;
        else                rfc_ctx.counting_method = RFC_COUNTING_METHOD_4PTM;
#elif RFC_HCM_SUPPORT
        if( use_hcm )       rfc_ctx.counting_method = RFC_COUNTING_METHOD_HCM;
        else                rfc_ctx.counting_method = RFC_COUNTING_METHOD_4PTM;
#elif RFC_ASTM_SUPPORT
        if( use_astm )      rfc_ctx.counting_method = RFC_COUNTING_METHOD_ASTM;
        else                rfc_ctx.counting_method = RFC_COUNTING_METHOD_4PTM;
#else /*!(RFC_HCM_SUPPORT  || RFC_ASTM_SUPPORT)*/
#if !RFC_MINIMAL
        rfc_ctx.counting_method = RFC_COUNTING_METHOD_4PTM;
#endif /*!RFC_MINIMAL*/
#endif /*(RFC_HCM_SUPPORT  || RFC_ASTM_SUPPORT)*/

        ok = RFC_feed( &rfc_ctx, buffer, data_len ) &&
             RFC_finalize( &rfc_ctx, residual_method );

        /* Free temporary buffer (cast) */
        if( (void*)buffer != (void*)data )
        {
            buffer = mem_alloc( buffer, 0, 0, RFC_MEM_AIM_TEMP );
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
                    mexErrMsgTxt( "Unexpected error occurred!" );
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
                        int i;
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
                            *ptr++ = (double)rfc_ctx.rfm[ MAT_OFFS( from, to ) ] / rfc_ctx.full_inc;
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
                        *ptr++ = (double)rfc_ctx.rp[i] / rfc_ctx.curr_inc;
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
                    /* assert( D == rfc_ctx.damage ); */
                    plhs[5] = tp;
                }
            }
#endif /*RFC_TP_SUPPORT*/
#if RFC_DH_SUPPORT
            /* Turning points */
            if( nlhs > 6 )
            {
                if( rfc_ctx.dh )
                {
                    mxArray *dh  = mxCreateDoubleMatrix( rfc_ctx.internal.pos, 1, mxREAL );
                    double  *dh_ptr = dh ? mxGetPr(dh) : NULL;

                    if( dh_ptr )
                    {
                        size_t i;

                        for( i = 0; i < rfc_ctx.internal.pos; i++ )
                        {
                            *dh_ptr++ = rfc_ctx.dh[i];
                        }
                    }

                    plhs[6] = dh;
                }
                else
                {
                    plhs[6] = mxCreateDoubleMatrix( 0, 0, mxREAL );
                }
            }
#endif /*RFC_DH_SUPPORT*/
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

        tp = (rfc_value_tuple_s*)CALLOC( mxGetNumberOfElements( mxData ), sizeof(rfc_value_tuple_s) );
        if( !tp )
        {
            mexErrMsgTxt( "Memory allocation error!" );
        }

        if( !RFC_init( &ctx, 0 /*class_count*/, 0.0 /*class_width*/, 0.0 /*class_offset*/, hysteresis, RFC_FLAGS_DEFAULT ) )
        {
            FREE(tp);
            mexErrMsgTxt( "Error on RFC init!" );
        }

        if( !RFC_tp_init( &ctx, tp, mxGetNumberOfElements( mxData ), true /*tp_is_static*/ ) )
        {
            FREE(tp);
            mexErrMsgTxt( "Error on RFC tp init!" );
        }

        if( !RFC_feed( &ctx, mxGetPr( mxData ), mxGetNumberOfElements( mxData ) ) )
        {
            FREE(tp);
            mexErrMsgTxt( "Error on RFC feed!" );
        }

        /*if( !finalize_res_ignore( &ctx, ctx.internal.flags ) )*/
        if( !RFC_finalize( &ctx, RFC_RES_IGNORE ) )
        {
            FREE(tp);
            mexErrMsgTxt( "Error on RFC finalize!" );
        }

        mxResult = mxCreateDoubleMatrix( 2, ctx.tp_cnt, mxREAL );
        mxAssert( mxResult, "Memory allocation error!" );
        for( ptr = mxGetPr( mxResult ), n = 0; n < ctx.tp_cnt; n++ )
        {
            *ptr++ = tp[n].value;
            *ptr++ = (double)tp[n].pos;
        }
        FREE( tp );

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
