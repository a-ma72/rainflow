
/* Using Rainflow C-Library in a C++ context */

#define RFC_TP_STORAGE std::vector<RF::rfc_value_tuple_s>

#include "../rainflow.hpp"
//#include "../neolib/segmented_array.hpp"

#define NUMEL(x) (sizeof(x)/sizeof(*(x)))


/* Made extern "C", so that C-Code (rfc_test.c) can access this function (in a C++ module)
   C++ functions are usually "name mangeled", that is averted here! */
extern "C"
bool wrapper_test_advanced( void )
{
    Rainflow rf;

    double data[] = { 1,6,2,8 };
    std::vector<double> data2;

    data2.push_back( 1 );
    data2.push_back( 6 );
    data2.push_back( 2 );
    data2.push_back( 8 );

#if 0
    rf.Parametrize( /*range_min*/ 0.0, /*range_max*/ 10.0,
                    /*range_fixed_min*/ 0.0, /*range_fixed_max*/ 10.0, /*range*/ DBL_NAN, 
                    /*class_count*/ 10, /*class_width*/ 1, 
                    /*dilation*/ 0.0, /*hysteresis*/ 0.5 );

    rf.DoRainflow( data, 4, true );
#endif

    rf.init( 100, 1, 0, 1 );
    rf.feed( data, NUMEL(data) );
    rf.feed( data2 );
    rf.finalize( Rainflow::RFC_RES_REPEATED );
    rf.deinit();

    return true;
}
