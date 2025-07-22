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


if (WIN32)
    set(PREBUILD_NUMPY_VERSIONS 1.23.2  2.0.0)  # Python 3.11
    # set(PREBUILD_NUMPY_VERSIONS 1.19.3  2.0.0)  # Python 3.9
else ()
    set(PREBUILD_NUMPY_VERSIONS)
endif ()
# Define NumPy versions used as prebuilds (used for MS Windows only), specifiers like ~=1.19.0 allowed.

# set (RFC_NUMPY_VERSION "1.23.5")
# NumPy version for the main module is defined per `RFC_NUMPY_VERSION`

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

    if (RFC_NUMPY_VERSION)
        numpy_pip_arg(PIP_ARG ${RFC_NUMPY_VERSION})
        execute_process(
                COMMAND "${Python3_EXECUTABLE}" -m pip install ${PIP_ARG}
        )
        # Install a specific version of NumPy if RFC_NUMPY_VERSION is defined.
    else ()
        execute_process(
                COMMAND "${Python3_EXECUTABLE}" -m pip install oldest-supported-numpy
        )
        # Otherwise, install the oldest supported version of NumPy.
    endif ()

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

# Ensure NumPy is installed, with the NumPy C API to build rfcnt
if (NOT Python3_NumPy_FOUND)
    message(STATUS "[Python] Didn't find NumPy, trying alternative search...")
    execute_process(
            COMMAND "${Python3_EXECUTABLE}" ${CMAKE_CURRENT_SOURCE_DIR}/cmake/numpy_get_include.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
            OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIRS
            RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )
    # Try to find NumPy include directories using a Python script if the initial search fails.

    if (Python3_NumPy_NOTFOUND)
        message(FATAL_ERROR "[Python] NumPy headers not found")
    else ()
        message(STATUS "[Python] NumPy found.")
    endif ()
endif ()

# Get the active NumPy version and directories
execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.__version__, end='')"
        OUTPUT_VARIABLE Python3_NumPy_VERSION
        RESULT_VARIABLE Python3_NumPy_NOTFOUND
)
# Retrieve the version of the installed NumPy package.

numpy_get_capi_version(Python3_NumPy_CAPI)

message(STATUS "[Python] Executable: ${Python3_EXECUTABLE}")
message(STATUS "[Python] Include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "[Python] Library: ${Python3_LIBRARIES}")
message(STATUS "[Python] NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "[Python] NumPy version: ${Python3_NumPy_VERSION}")
message(STATUS "[Python] NumPy API version: ${Python3_NumPy_CAPI}")
# Output status messages with details about the Python3 and NumPy setup.


# ======================================== Collect NumPy include directories ==========================================
# Unfortunately, CMake does not provide a way to change the compilers include directories at build time.
# Therefore all packages must be installed at configure time in order to find out the directories to pass them
# to compiler options.
message(STATUS "Collect NumPy include directories")
# Output a status message indicating the start of collecting NumPy include directories.

file(REMOVE "prebuilds.json")
if (EXISTS "prebuilds.json")
    message(FATAL_ERROR "Could not delete 'prebuilds.json'")
endif ()
set(PREBUILDS_JSON "")
set(JSON_SEP " ")

set(NUMPY_VERSIONS_LIST ${PREBUILD_NUMPY_VERSIONS})
# Collect numpy versions needed into NUMPY_VERSIONS_LIST

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/_ext)
# Ensure directory for prebuilds

# Initialize the list of NumPy versions with the given PREBUILD_NUMPY_VERSIONS.
set(NUMPY_INCLUDE_DIRECTORIES)
# Initialize an empty list to store the NumPy include directories.

list(APPEND NUMPY_VERSIONS_LIST oldest-supported-numpy ${RFC_NUMPY_VERSION})
# Append the "oldest-supported-numpy" and RFC_NUMPY_VERSION to the list of NumPy versions.

foreach (NUMPY_VERSION ${NUMPY_VERSIONS_LIST})
    # Loop through each NumPy version in the list.

    numpy_pip_arg(PIP_ARG ${NUMPY_VERSION})

    execute_process(
            COMMAND "${Python3_EXECUTABLE}" -mpip uninstall -y oldest-supported-numpy numpy
    )
    # Uninstall any existing NumPy versions to ensure a clean installation.

    execute_process(
            COMMAND "${Python3_EXECUTABLE}" -mpip install ${PIP_ARG}
            COMMAND_ECHO STDOUT
    )
    # Install the specified version of NumPy.

    execute_process(
            COMMAND "${Python3_EXECUTABLE}" ${CMAKE_CURRENT_SOURCE_DIR}/cmake/numpy_get_include.py
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
            OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIR
            RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )
    # Run a Python script to get the include directory for NumPy and store it in Python3_NumPy_INCLUDE_DIR.

    if (Python3_NumPy_NOTFOUND)
        message(FATAL_ERROR "NumPy not found")
        # If the include directory was not found, output an error message and terminate the process.
    endif ()

    # Get the active NumPy version and directories
    execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.__version__, end='')"
            OUTPUT_VARIABLE Python3_NumPy_VERSION
    )
    # Retrieve the version of the installed NumPy package.

    numpy_get_capi_version(NUMPY_CAPI_VERSION)
    prebuild_target_name(PREBUILD_TARGET ${Python3_NumPy_VERSION} ${NUMPY_CAPI_VERSION})
    string(JSON VERSION_EXISTS ERROR_VARIABLE ERROR_MSG GET "{${PREBUILDS_JSON}}" "${Python3_NumPy_VERSION}")
    if (ERROR_MSG)
        string(APPEND PREBUILDS_JSON "    ${JSON_SEP}\"${Python3_NumPy_VERSION}\": {\"capi\": ${NUMPY_CAPI_VERSION}, \"target\": \"${PREBUILD_TARGET}\"}\n")
        set(JSON_SEP ",")
        list(APPEND NUMPY_INCLUDE_DIRECTORIES  ${Python3_NumPy_INCLUDE_DIR})
    endif ()
    # Append the found include directory to the list of NumPy include directories.
endforeach ()

string(APPEND PREBUILDS_JSON "    ${JSON_SEP}\"${Python3_NumPy_VERSION}\": {\"capi\":${NUMPY_CAPI_VERSION}, \"target\": \"root\"}\n")
file(WRITE "prebuilds.json" "{\n${PREBUILDS_JSON}}")
file(READ "prebuilds.json" PREBUILDS_JSON)
list(REMOVE_DUPLICATES NUMPY_INCLUDE_DIRECTORIES)
# Remove duplicate entries from the list of NumPy include directories.

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.c"
            "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/../lib/rainflow.hpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/../lib/config.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)
# Copy rainflow core sources


# ========================================= Create main module and prebuilds ==========================================

# Add main module with the current installed numpy version.
rfcnt_library("${Python3_NumPy_VERSION}" ""  "${PROJECT_NAME}"  "."  rfcnt_module)

# Add prebuilds for different NumPy packages
set(RFCNT_PREBUILD_TARGETS rfcnt_module)
# Initialize an empty list to store prebuild targets.
# The list is used as a target dependency for consecutive targets to ensure build order and
# prevent conflicts when concurrent processes would update the virtual environment.

string(JSON PREBUILDS_COUNT LENGTH ${PREBUILDS_JSON})
foreach (nr RANGE 1 ${PREBUILDS_COUNT})
    # Loop through each version of NumPy in the list of PREBUILD_NUMPY_VERSIONS.

    math (EXPR member_id "${nr} - 1")
    string(JSON VERSION_STRING MEMBER ${PREBUILDS_JSON} ${member_id})
    string(JSON PREBUILD_NAME GET ${PREBUILDS_JSON} ${VERSION_STRING} "target")

    if (NOT ${PREBUILD_NAME} STREQUAL "root")
        message(STATUS ${PREBUILD_NAME})
        # Output a status message with the prebuild name.

        set(RFCNT_PREBUILD_TARGET target_${PREBUILD_NAME})
        # Define the prebuild target name based on the prebuild name.

        if (NOT TARGET ${RFCNT_PREBUILD_TARGET})
            # Check if the target does not already exist.

            rfcnt_library(${VERSION_STRING} "${RFCNT_PREBUILD_TARGETS}"  ${PREBUILD_NAME}  _ext  ${RFCNT_PREBUILD_TARGET})
            # Call the rfcnt_library function to create the library target with the specified details.

            list(APPEND RFCNT_PREBUILD_TARGETS ${RFCNT_PREBUILD_TARGET})
            # Append the new prebuild target to the list of prebuild targets.
        endif ()
    endif ()
endforeach ()


add_custom_target(
        rfcnt_prebuilds
        DEPENDS ${RFCNT_PREBUILD_TARGETS}
        COMMENT "Build rfcnt prebuilds"
)
# Custom target for prebuilds, preserve the build order.


add_custom_target(
        rfcnt
        DEPENDS rfcnt_prebuilds rfcnt_module
        COMMENT "Build rfcnt modules"
)
# Custom target for main module and prebuilds.


add_custom_target(
        build_wheel_isolated
        COMMAND "${CMAKE_COMMAND}" -E remove -f $<TARGET_FILE:rfcnt_module>
        COMMAND "${Python3_EXECUTABLE}" -mbuild
        DEPENDS rfcnt
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Build rfcnt wheel (needs a Numpy package defined in pyproject.toml!)"
)
# Custom target for creating sdist and wheel from sdist for rfcnt with setuptools (prebuilds by CMake).


add_custom_target(
        build_wheel
        COMMAND "${Python3_EXECUTABLE}" -m pip uninstall -y numpy
        COMMAND "${Python3_EXECUTABLE}" -m pip install oldest-supported-numpy
        COMMAND "${Python3_EXECUTABLE}" -mbuild --no-isolation
        DEPENDS rfcnt
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        COMMENT "Build rfcnt wheel"
)
# Custom target for creating wheel with main module and prebuilds created by CMake.


add_custom_target(
        clean_all
        COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/build"
        COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/dist"
        COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/rfcnt.egg-info"
        COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/venv"
        COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_CURRENT_SOURCE_DIR}/_ext"
        COMMAND "${CMAKE_COMMAND}" -E remove "${CMAKE_CURRENT_SOURCE_DIR}/prebuilds.json"
)
# Custom target for cleanup directories.
