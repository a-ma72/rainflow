# Function wrapper for MATLAB(R)

# Regular build:
# ===============
#[[
    cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_EXPORT_MEX=1
    cmake --build build/src/matlab --target rfc_mex --config Release
#]]

# Build using special MATLAB(R) installation:
# ===========================================
#[[
    cmake -S. -Bbuild -DRFC_EXPORT_PY=0 -DRFC_EXPORT_MEX=1 -DMatlab_ROOT_DIR=/appl/matlab/2017b
    cmake --build build/src/matlab --target rfc_mex --config Release
#]]

# Build packages
# ==============
#[[
    cd build
    cpack --config CPackConfig.cmake
    cpack --config CPackSourceConfig.cmake
#]]


cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0048 NEW)  # The project() command manages VERSION variables.
if (POLICY CMP0076)
    cmake_policy(SET CMP0076 NEW)  # The target_sources() command converts relative paths to absolute.
endif ()
cmake_policy(SET CMP0091 NEW)  # MSVC runtime library flags are selected by an abstraction.

# Project name
project(rfc_mex LANGUAGES C CXX)

# find MATLAB
# (https://de.mathworks.com/support/requirements/previous-releases.html)
set(MATLAB_FIND_DEBUG 1)
set(MATLAB_ADDITIONAL_VERSIONS
        "R2019b=9.7"
        "R2020a=9.8"
        "R2020b=9.9"
        "R2021a=9.10"
        "R2021b=9.11"
        "R2022a=9.12"
        "R2022b=9.13"
        "R2023a=9.14"
        "R2023b=9.15"
        "R2024a=9.16"
        "R2024b=9.17"
)
set(VERSION_MATLAB "" CACHE STRING "The desired Matlab installation (`9.7` for example), or autodetect if empty.")
if (DEFINED ENV{Matlab_ROOT_DIR})
    set(Matlab_ROOT_DIR $ENV{Matlab_ROOT_DIR})
endif ()
if ("${VERSION_MATLAB}" STREQUAL "")
    find_package(Matlab COMPONENTS MAIN_PROGRAM)
else ()
    find_package(Matlab ${VERSION_MATLAB} EXACT COMPONENTS MAIN_PROGRAM)
endif ()

if (NOT MATLAB_FOUND)
    message(WARNING "MATLAB not found, unable to build MEX file.")
else ()
    message(STATUS "MATLAB found, able to build MEX file.")
    if (UNIX)
        matlab_get_version_from_matlab_run(${Matlab_MAIN_PROGRAM} matlab_version)
    else ()
        set(matlab_version ${VERSION_MATLAB})
    endif ()

    # MEX function (MATLAB)
    if (RFC_EXPORT_MEX)
        matlab_add_mex(
                NAME rfc_mex
                SRC ${CMAKE_CURRENT_SOURCE_DIR}/rfc.c
                LINK_TO ${LIBM_LIBRARY}
                rfc_core
                OUTPUT_NAME rfc
        )
    else ()
        matlab_add_mex(
                NAME rfc_mex
                SRC ${CMAKE_CURRENT_SOURCE_DIR}/rfc.c
                ${rfc_core_SOURCE_DIR}/rainflow.c
                LINK_TO ${LIBM_LIBRARY}
                OUTPUT_NAME rfc
        )
        target_include_directories(rfc_mex PRIVATE
                ${rfc_core_SOURCE_DIR}
        )
        target_compile_definitions(rfc_mex PRIVATE
                RFC_VALUE_TYPE=double
                # RFC_USE_INTEGRAL_COUNTS
                RFC_USE_DELEGATES
                RFC_USE_HYSTERESIS_FILTER
                RFC_GLOBAL_EXTREMA
                RFC_HCM_SUPPORT
                RFC_ASTM_SUPPORT
                RFC_TP_SUPPORT
                RFC_DH_SUPPORT
                RFC_AT_SUPPORT
                RFC_AR_SUPPORT
                RFC_DAMAGE_FAST
                # RFC_DEBUG_FLAGS
                RFC_EXPORT_MEX
                # RFC_UNIT_TEST
        )
    endif ()

    # install to source root by default
    install(TARGETS rfc_mex LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
endif ()


# CPack
include(InstallRequiredSystemLibraries)

set(CPACK_GENERATOR TGZ ZIP)
set(CPACK_SOURCE_GENERATOR TGZ ZIP)
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/package)
#set( CPACK_PACKAGE_NAME "rainflow")
set(CPACK_PACKAGE_VERSION_MAJOR ${RFC_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${RFC_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH "")
set(CPACK_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../../README.md")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fast rainflow counting written in C (C99)")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE")
set(CPACK_SOURCE_IGNORE_FILES
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
set(CPACK_STRIP_FILES TRUE)
set(CPACK_SOURCE_STRIP_FILES TRUE)

include(CPack)
