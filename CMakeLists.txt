# First set Matlab_ROOT_DIR environment variable to your installed matlab path,
# such as 'export Matlab_ROOT_DIR=/usr/local/MATLAB/R2017b'!

# Building makefiles:
# mkdir build
# cd build
# cmake -G "Visual Studio 14 2015 Win64" ..


cmake_minimum_required( VERSION 3.2 )
set(CMAKE_CONFIGURATION_TYPES Release Debug )
project( Rainflow )

set( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/build )

set( MATLAB_FIND_DEBUG 1 )
find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY )

IF( MATLAB_FOUND )
    message( STATUS "MATLAB Found, MATLAB MEX will be compiled." )
    message( STATUS ${Matlab_LIBRARIES} )
ELSE( MATLAB_FOUND )
    MESSAGE( "MATLAB not found...nothing will be built." )
ENDIF( MATLAB_FOUND )

add_definitions( -DMATLAB_MEX_FILE ) #define matlab macros

# set up matlab libraries
INCLUDE_DIRECTORIES( ${Matlab_INCLUDE_DIRS} )

matlab_add_mex( NAME Rainflow SRC rainflow.c OUTPUT_NAME rfc )
target_link_libraries( Rainflow ${Matlab_LIBRARIES} )

SET_TARGET_PROPERTIES( Rainflow PROPERTIES SUFFIX .${Matlab_MEX_EXTENSION})

# install to /mex by default
add_custom_command( TARGET Rainflow POST_BUILD 
                    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/mex/$<CONFIG>/
                    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Rainflow> ${CMAKE_INSTALL_PREFIX}/mex/$<CONFIG>/ )
