# First set Matlab_ROOT_DIR environment variable to your installed matlab path,
# such as 'export Matlab_ROOT_DIR=/usr/local/MATLAB/R2017b'!
# (GNUC: gcc for C file, g++ for CXX files)

# Building makefiles:
# mkdir build
# cd build
# cmake -G "Visual Studio 14 2015 Win64" ..


cmake_minimum_required( VERSION 3.2 )
set( CMAKE_CONFIGURATION_TYPES Release Debug )
project( Rainflow )
set( RFC_VERSION_MAJOR "0" )
set( RFC_VERSION_MINOR "1" )

option( RFC_USE_INTEGRAL_COUNTS "Use integral (not floating) data type for counts" OFF )
option( RFC_USE_DELEGATES       "Use delegates (functors)" OFF )
option( RFC_GLOBAL_EXTREMA      "Always calculate global extrema" OFF )
option( RFC_HCM_SUPPORT         "Support HCM (Clormann/Seeger) algorithm" ON )
option( RFC_TP_SUPPORT          "Support turning points" ON )


set( CMAKE_POSITION_INDEPENDENT_CODE ON )

if( MSVC )
    # Turn off misguided "secure CRT" warnings in MSVC.
    # Microsoft wants people to use the MS-specific <function>_s
    # versions of certain C functions but this is difficult to do
    # in platform-independent code.
    add_definitions( -D_CRT_SECURE_NO_WARNINGS )
endif( MSVC )

set( MATLAB_FIND_DEBUG 1 )
find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY )

if( MATLAB_FOUND )
    message( STATUS "MATLAB Found, MATLAB MEX will be compiled." )
    # message( STATUS ${Matlab_LIBRARIES} )
else()
    message( FATAL_ERROR "MATLAB not found...nothing will be built." )
endif()

add_definitions( -DMATLAB_MEX_FILE )

# set up matlab libraries
include_directories( ${Matlab_INCLUDE_DIRS} )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/config.h )

matlab_add_mex( NAME ${CMAKE_PROJECT_NAME} SRC rainflow.c OUTPUT_NAME rfc )
target_link_libraries( ${CMAKE_PROJECT_NAME} ${Matlab_LIBRARIES} )

# install to /bin and /share by default
install( TARGETS ${CMAKE_PROJECT_NAME} RUNTIME DESTINATION bin LIBRARY DESTINATION bin )


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
     "/\\\\.build/"

     # Temporary files.
     "\\\\.#"
     "/#"
     "~$"
    )
set( CPACK_STRIP_FILES TRUE )
set( CPACK_SOURCE_STRIP_FILES TRUE )

include( CPack )
