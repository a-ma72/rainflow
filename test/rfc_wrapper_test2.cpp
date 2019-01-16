
/* Using Rainflow C-Library in a C++ context */

#include "../rainflow.hpp"


/* Made extern "C", so that C-Code (rfc_test.c) can access this function (in a C++ module)
   C++ functions are usually "name mangeled", that is averted here! */
extern "C"
bool wrapper_test_advanced( void )
{
    CRainflow rf;
    double data[] = { 1,6,2,8 };

    rf.Parametrize( /*range_min*/ 0.0, /*range_max*/ 10.0,
                    /*range_fixed_min*/ 0.0, /*range_fixed_max*/ 10.0, /*range*/ DBL_NAN, 
                    /*class_count*/ 10, /*class_width*/ 1, 
                    /*dilation*/ 0.0, /*hysteresis*/ 0.5 );

    rf.DoRainflow( data, 4, true );

    return true;
}
