# CMakeLists for project "rainflow"
# 2023, Andreas Martin

# See `README.md` for build and installation hints.


cmake_minimum_required(VERSION 3.14)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_BUILD_PARALLEL_LEVEL 12)

if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()

# Project name and version
project(rainflow VERSION 4.7 LANGUAGES NONE)

# C++11
enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug) 
    message(STATUS "Rainflow: Debug build") 
    add_definitions(-DDEBUG -D_DEBUG) 
endif()

# Compiler dependencies
if(MSVC)
    # Turn off misguided "secure CRT" warnings in MSVC.
    # Microsoft wants people to use the MS-specific <function>_s
    # versions of certain C functions but this is difficult to do
    # in platform-independent code.
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Math library
find_library(LIBM_LIBRARY NAMES m)
if(NOT LIBM_LIBRARY)
    set(LIBM_LIBRARY "")
endif()

# Options valid, if project compiled as standalone only
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(STATUS "Build ${PROJECT_NAME} as main project")
    option(RFC_MINIMAL                "Build module with as few features as possible" OFF)
    option(RFC_USE_INTEGRAL_COUNTS    "Use integral (non-floating) data type for counts" OFF)
    option(RFC_USE_HYSTERESIS_FILTER  "Use hysteresis filtering" ON)
    option(RFC_USE_DELEGATES          "Use delegates (functors)" ON)
    option(RFC_GLOBAL_EXTREMA         "Always calculate global extrema" ON)
    option(RFC_HCM_SUPPORT            "Support HCM (Clormann/Seeger) algorithm" ON)
    option(RFC_ASTM_SUPPORT           "Support ASTM E 1049-85 algorithm" ON)
    option(RFC_TP_SUPPORT             "Support turning points" ON)
    option(RFC_DH_SUPPORT             "Support \"spread damage\" over turning points and damage history" ON)
    option(RFC_AT_SUPPORT             "Support amplitude transformation regarding mean load influence on fatigue strength" ON)
    option(RFC_AR_SUPPORT             "Support automatic growth of counting buffers" ON)
    option(RFC_DAMAGE_FAST            "Enables fast damage calculation (per look-up table)" ON)
    option(RFC_DEBUG_FLAGS            "Enables flags for detailed examination" OFF)
    option(RFC_EXPORT_MEX             "Export a function wrapper for MATLAB(R)" ON)
    option(RFC_EXPORT_PY              "Export a function wrapper for Python)" ON)
    option(RFC_UNIT_TEST              "Generate rainflow testing program for unit test" ON)
    set(RFC_VALUE_TYPE      double    CACHE STRING "Value type of input data to be processed")
    set(RFC_PYTHON_VERSION  "3.9"     CACHE STRING "Expected Python version")
    set(RFC_NUMPY_VERSION   ""        CACHE STRING "Numpy version to link to")
else()
    message(STATUS "Build ${PROJECT_NAME} as subsequent project")
endif()

set(RFC_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(RFC_VERSION_MINOR ${PROJECT_VERSION_MINOR})

# Save options in configuration file
configure_file(${CMAKE_CURRENT_LIST_DIR}/src/lib/config.h.in ${CMAKE_CURRENT_LIST_DIR}/src/lib/config.h)

# Targets
add_subdirectory(src)

# Unit tests
if(${RFC_UNIT_TEST})
    enable_testing()
    include(CTest)
    add_subdirectory(test)  # EXCLUDE_FROM_ALL)
    add_test(NAME rfc_unit_test COMMAND rfc_test)
endif()
