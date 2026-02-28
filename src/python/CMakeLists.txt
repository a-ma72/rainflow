# Function wrapper for Python

# Regular build
# =============
#[[
    cmake -S. -Bbuild -DRFC_EXPORT_PY=1 -DRFC_EXPORT_MEX=0
    cmake --build build --target build_wheel_isolated --config Release

    # MinGW
    cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ ...
#]]

# Linux build
# ===========
#[[
    rm -rf build
    cmake -G "Unix Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ ...
    cmake --build build --target build_wheel_isolated --config Release
#]]

# Python Setup
#[[
    pip install --force-reinstall --no-deps package.tar
    # mingw32:
    # in PYTHONPATH\Lib\distutils, create a file distutils.cfg and add these lines:
    [build]
    compiler = mingw32
    # Or:
    pip install --global-option build_ext --global-option --compiler=mingw32 <package_zip>
#]]



cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0048 NEW)  # The project() command manages VERSION variables.
if (POLICY CMP0076)
    cmake_policy(SET CMP0076 NEW)  # The target_sources() command converts relative paths to absolute.
endif ()
cmake_policy(SET CMP0091 NEW)  # MSVC runtime library flags are selected by an abstraction.
cmake_policy(SET CMP0094 NEW)  # Modules FindPython3, FindPython2 and FindPython use LOCATION for lookup strategy.

set (Python3_FIND_REGISTRY "LAST")
# Search the Windows registry for Python installations last

set (Python3_FIND_VIRTUALENV "FIRST")
# Prefer to find Python in a virtual environment if available

# set (CMAKE_FIND_DEBUG_MODE TRUE)
# Uncommenting this line enables debug mode for CMake, which provides detailed information about the find process.

project(rfcnt LANGUAGES C CXX)
# Project name

if (UNIX)
    set(PATH_SEP ":")
    set(PYTHON_PLATFORM "linux_x86_64")
else ()
    set(PATH_SEP "\;")
    set(PYTHON_PLATFORM "win_amd64")
endif ()

include (${CMAKE_CURRENT_SOURCE_DIR}/cmake/rfcnt_target_functions.cmake)


# =============================================== Find Python 3 & NumPy ===============================================
set(Python3_FIND_ABI "OFF" "ANY" "ANY")
# This setting ensures that CMake does not consider the ABI when finding Python3, making the search less restrictive.

if (VIRTUAL_ENV)
    # Prefer active virtual environment
    message(STATUS "Using virtual environment ${VIRTUAL_ENV} (cmake argument preset)")
    set(Python3_ROOT_DIR ${VIRTUAL_ENV})
elseif (DEFINED ENV{VIRTUAL_ENV})
    # Prefer active virtual environment
    message(STATUS "Using active virtual environment $ENV{VIRTUAL_ENV}")
    set(Python3_ROOT_DIR $ENV{VIRTUAL_ENV})
else ()
    message(STATUS "(Re-)create virtual environment")
    # If no virtual environment is active or specified, create a new one.

    file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/venv")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/venv")
        message(FATAL_ERROR "Could not delete local virtual environment '${CMAKE_CURRENT_SOURCE_DIR}/venv'")
    endif ()
    # Remove any existing virtual environment in the source directory.

    find_package(Python3 ${RFC_PYTHON_VERSION} REQUIRED COMPONENTS Interpreter)
    # Find the Python3 interpreter first.

    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m venv ${CMAKE_CURRENT_SOURCE_DIR}/venv
    )
    # Create a new virtual environment in the source directory.

    set(Python3_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/venv")
    # Set the root directory for Python3 to the new virtual environment.
    if (WIN32)
        set(Python3_EXECUTABLE "${Python3_ROOT_DIR}/Scripts/python.exe")
    else ()
        set(Python3_EXECUTABLE "${Python3_ROOT_DIR}/bin/python3")
    endif ()
    # Preliminary redirect to the just created environment

    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m pip install -U pip setuptools wheel build
    )
    # Install necessary Python packages in the virtual environment.

    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m pip install "numpy>=2.0"
    )
    # Install NumPy 2.x for building with NPY_TARGET_VERSION backward compatibility.

    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
        OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIR
        RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )
    # Get the include directory for NumPy.

    unset(Python3_INCLUDE_DIRS CACHE)
    unset(Python3_LIBRARIES CACHE)
    unset(Python3_NumPy_INCLUDE_DIRS CACHE)
    unset(Python3_NumPy_VERSION CACHE)
    # Unset variables to ensure a clean slate for the next find_package call.
endif ()
find_package(Python3 ${RFC_PYTHON_VERSION} COMPONENTS Development Interpreter Numpy)
# Find the Python3 package with the specified version and required components.

# Fallback: get Python include dirs if find_package failed
if (NOT Python3_INCLUDE_DIRS)
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('include'), end='')"
        OUTPUT_VARIABLE Python3_INCLUDE_DIRS
    )
    string(STRIP "${Python3_INCLUDE_DIRS}" Python3_INCLUDE_DIRS)
    message(STATUS "[Python] Include dirs (fallback): ${Python3_INCLUDE_DIRS}")
endif ()

# Ensure NumPy is installed, with the NumPy C API to build rfcnt
if (NOT Python3_NumPy_FOUND OR NOT Python3_NumPy_INCLUDE_DIRS)
    message(STATUS "[Python] Didn't find NumPy, trying alternative search...")
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
        OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIRS
        RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )
    # Try to find NumPy include directories using a direct Python command if find_package fails.

    if (Python3_NumPy_NOTFOUND)
        message(FATAL_ERROR "[Python] NumPy headers not found")
    else ()
        string(STRIP "${Python3_NumPy_INCLUDE_DIRS}" Python3_NumPy_INCLUDE_DIRS)
        message(STATUS "[Python] NumPy found (fallback): ${Python3_NumPy_INCLUDE_DIRS}")
    endif ()
endif ()

# Get the active NumPy version and directories
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.__version__, end='')"
    OUTPUT_VARIABLE Python3_NumPy_VERSION
    RESULT_VARIABLE Python3_NumPy_NOTFOUND
)
# Retrieve the version of the installed NumPy package.

message(STATUS "[Python] Executable: ${Python3_EXECUTABLE}")
message(STATUS "[Python] Include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "[Python] Library: ${Python3_LIBRARIES}")
message(STATUS "[Python] NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "[Python] NumPy version: ${Python3_NumPy_VERSION}")
# Output status messages with details about the Python3 and NumPy setup.


# ============================================= Copy core sources and build ============================================

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/../lib/config.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)
# Copy rainflow core sources

# Set rfc_core_SOURCE_DIR so rfcnt_target_functions.cmake can find the headers
set(rfc_core_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")


# ========================================= Create main module =========================================================

# Add main module with the current installed numpy version.
rfcnt_library("${PROJECT_NAME}"  "."  rfcnt_module)


add_custom_target(
    rfcnt
    DEPENDS rfcnt_module
    COMMENT "Build rfcnt module"
)
# Custom target for main module.


if(CMAKE_GENERATOR MATCHES "MinGW Makefiles")
    set(RFCNT_BUILD_COMPILER_SETTING --config-setting="--build-option=build_ext" --config-setting="--build-option=--compiler=mingw32")
else()
    set(RFCNT_BUILD_COMPILER_SETTING "")
endif()

add_custom_target(
    build_wheel_isolated
    COMMAND "${CMAKE_COMMAND}" -E remove -f $<TARGET_FILE:rfcnt_module>
    COMMAND "${Python3_EXECUTABLE}" -m build --no-isolation --wheel
        ${RFCNT_BUILD_COMPILER_SETTING}
        --outdir "${CMAKE_CURRENT_SOURCE_DIR}/dist"
    DEPENDS rfcnt_module
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Build rfcnt wheel (PEP 517)"
)
# Custom target for creating wheel via PEP 517 build frontend.


add_custom_target(
    build_wheel
    COMMAND "${Python3_EXECUTABLE}" -mbuild --no-isolation
    DEPENDS rfcnt_module
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Build rfcnt wheel"
)
# Custom target for creating wheel with main module.


add_custom_target(
    clean_all
    COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/build"
    COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/dist"
    COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/rfcnt.egg-info"
    COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/venv"
)
# Custom target for cleanup directories.
