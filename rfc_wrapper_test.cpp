
/* Using Rainflow C-Library in an c++ context */

namespace rainflow_C
{
    extern "C"
    {
        #include "rainflow.h"
    }
}

extern "C"
bool wrapper_test( void )
{
    namespace rf = rainflow_C;

    rf::rfc_ctx ctx = { sizeof(ctx) };

    return
    rf::RFC_init( &ctx, /*class_count*/  10, 
                        /*class_width*/  1, 
                        /*class_offset*/ 0,
                        /*hysteresis*/   1,
                        /*flags*/        rf::RFC_FLAGS_DEFAULT )
    &&
    rf::RFC_deinit( &ctx );
}
