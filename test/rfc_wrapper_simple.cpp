
/* Using Rainflow C-Library in a C++ context */

#include "../rainflow.h"


/* Made extern "C", so that C-Code (rfc_test.c) can access this function (in a C++ module)
   C++ functions are usually "name mangeled", that is averted here! */
extern "C"
bool wrapper_test_simple( void )
{
    namespace rf = rainflow_C;  /* Using a namespace alias to make it short */

    rf::rfc_ctx_s ctx = { sizeof(ctx) };

    /* Just init and deinit here */
    return
    rf::RFC_init( &ctx, /*class_count*/  10, 
                        /*class_width*/  1, 
                        /*class_offset*/ 0,
                        /*hysteresis*/   1,
                        /*flags*/        rf::RFC_FLAGS_DEFAULT )
    &&
    rf::RFC_deinit( &ctx );
}
