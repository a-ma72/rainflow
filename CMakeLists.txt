# First set Matlab_ROOT_DIR environment variable to your installed matlab path,
# such as 'export Matlab_ROOT_DIR=/usr/local/MATLAB/R2017b'!
# (GNUC: gcc for C file, g++ for CXX files)

# Building makefiles:
# mkdir build
# cd build
# cmake -G "Visual Studio 14 2015 Win64" ..


cmake_minimum_required( VERSION 3.6 )
project( Rainflow )
set( CMAKE_CONFIGURATION_TYPES Release Debug )
set( RFC_VERSION_MAJOR "0" )
set( RFC_VERSION_MINOR "1" )

# options valid, if project compiled as standalone only
if( CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( STATUS "Build ${PROJECT_NAME} as main project")
    option( RFC_USE_INTEGRAL_COUNTS "Use integral (non-floating) data type for counts" OFF )
else()
    message( STATUS "Build ${PROJECT_NAME} as subsequent project")
endif()

set( CMAKE_POSITION_INDEPENDENT_CODE ON )

if( MSVC )
    # Turn off misguided "secure CRT" warnings in MSVC.
    # Microsoft wants people to use the MS-specific <function>_s
    # versions of certain C functions but this is difficult to do
    # in platform-independent code.
    add_definitions( -D_CRT_SECURE_NO_WARNINGS )
endif()

if( UNIX )
    set( MATH_LIB m )
endif()

set( MATLAB_FIND_DEBUG 1 )
find_package( Matlab QUIET REQUIRED COMPONENTS MX_LIBRARY )

if( MATLAB_FOUND )
    message( STATUS "MATLAB Found, MATLAB MEX will be compiled." )
    # message( STATUS ${Matlab_LIBRARIES} )
else()
    message( FATAL_ERROR "MATLAB not found...nothing will be built." )
endif()

# set up matlab libraries
include_directories( ${Matlab_INCLUDE_DIRS} )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h )

matlab_add_mex( NAME ${PROJECT_NAME} SRC rainflow.c OUTPUT_NAME rfc )
target_compile_definitions( ${PROJECT_NAME} PRIVATE MATLAB_MEX_FILE _SCL_SECURE_NO_WARNINGS )
target_link_libraries( ${PROJECT_NAME} ${Matlab_LIBRARIES} ${MATH_LIB} )

add_library( rfc STATIC rainflow.c )
target_link_libraries( ${PROJECT_NAME} ${MATH_LIB} )

# Test application, start project for MSVC
add_executable( rfc_test rainflow.c test/rfc_test.c test/rfc_wrapper_simple.cpp )
target_compile_definitions( rfc_test PRIVATE _SCL_SECURE_NO_WARNINGS )
target_link_libraries( rfc_test ${MATH_LIB} )
target_sources( rfc_test PUBLIC greatest/greatest.h )
set_target_properties( rfc_test PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS YES )
set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT rfc_test )

# install to /bin by default
install( TARGETS ${PROJECT_NAME} rfc_test RUNTIME DESTINATION bin LIBRARY DESTINATION bin )
install( FILES test/long_series.csv test/long_series.c matlab/validate.m DESTINATION bin )


# CPack
include( InstallRequiredSystemLibraries )

set( CPACK_GENERATOR TGZ ZIP )
set( CPACK_SOURCE_GENERATOR TGZ ZIP )
set( CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/package )
set( CPACK_PACKAGE_VERSION_MAJOR ${RFC_VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${RFC_VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH "" )
set( CPACK_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR} )
set( CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.MD" )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast rainflow counting written in C (C99)" )
set( CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" )
set( CPACK_SOURCE_IGNORE_FILES 
     # Files specific to version control.
     "/\\\\.git/"
     "/\\\\.gitattributes$"
     "/\\\\.github/"
     "/\\\\.gitignore$"
     "/\\\\.hooks-config$"

     # Package build.
     "/build"
     "/.git*"

     # Temporary files.
     "\\\\.#"
     "/#"
     "~$"
    )
set( CPACK_STRIP_FILES TRUE )
set( CPACK_SOURCE_STRIP_FILES TRUE )

include( CPack )
