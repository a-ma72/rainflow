
/* Using Rainflow C-Library in a C++ context */

#include "../rainflow.hpp"


/* Made extern "C", so that C-Code (rfc_test.c) can access this function (in a C++ module)
   C++ functions are usually "name mangeled", that is averted here! */
extern "C"
bool wrapper_test2( void )
{
    CRainflow rf;

    rf.CalcDamage( 10, 10, 1 );

    return true;
}
