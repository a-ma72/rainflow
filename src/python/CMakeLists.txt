# Function wrapper for Python

# Regular build
# =============
#[[
    cmake -S. -Bbuild -DRFC_EXPORT_PY=1 -DRFC_EXPORT_MEX=0
    cmake --build build/src/python --target rfcnt --config Release

    # MinGW
    cmake -G "MinGW Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ ...
#]]

# Linux build
# ===========
#[[
    rm -rf build
    cmake -G "Unix Makefiles" -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ ...
    cmake --build build/src/python --target rfcnt --config Release
#]]

# Build Python Setup
#[[
    python -m build -nwx
    python setup.py bdist_wheel --plat-name=win-amd64
    python setup.py bdist --formats=wininst
    python setup.py bdist_wininst --title= --bitmap=
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

# Python Package for PyPi
#[[
    # Update PYTHONHOME environment variable
    # Update PATH environment variable
    # Make cmake settings (cmake-gui)
    # Open terminal in project root path:
        cmake -S. -Bbuild -DRFC_EXPORT_MEX=0 -DRFC_EXPORT_PY=1 -DRFC_PYTHON_VERSION="3.9"
        cd build
        py -mvenv venv
        venv\Scripts\activate.bat
        pip install setuptools wheel build oldest-supported-numpy
        cd ..\src\python
        make_prebuilds.bat  # (windows only)
        py -mbuild -nw
        cd ..\..
        build\venv\Scripts\deactivate.bat
#]]


cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0048 NEW)  # The project() command manages VERSION variables.
if (POLICY CMP0076)
    cmake_policy(SET CMP0076 NEW)  # The target_sources() command converts relative paths to absolute.
endif ()
cmake_policy(SET CMP0091 NEW)  # MSVC runtime library flags are selected by an abstraction.

# Project name
project(rfcnt LANGUAGES C CXX)
# set (RFC_NUMPY_VERSION "1.23.5")

if (UNIX)
    set(PATH_SEP ":")
    set(PYTHON_PLATFORM "linux_x86_64")
else ()
    set(PATH_SEP "\;")
    set(PYTHON_PLATFORM "win_amd64")
endif ()

# Find Python 3
set(Python3_FIND_ABI "OFF" "ANY" "ANY")
# set (CMAKE_FIND_DEBUG_MODE TRUE)
if (DEFINED ENV{VIRTUAL_ENV})
    # Prefer active virtual environment
    set(Python3_ROOT_DIR $ENV{VIRTUAL_ENV})
elseif (DEFINED ENV{PYTHONHOME})
    # Prefer predefined Python home path
    set(Python3_ROOT_DIR $ENV{PYTHONHOME})
endif ()
find_package(Python3 ${RFC_PYTHON_VERSION} COMPONENTS Development Interpreter Numpy)
if (NOT ("${RFC_NUMPY_VERSION}" STREQUAL ""))
    # Build module for pinned numpy version
    set(NUMPY_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/numpy/${RFC_NUMPY_VERSION})
    if (NOT EXISTS ${NUMPY_PACKAGE})
        # Local numpy package doesn't exist yet, install it
        execute_process(
                COMMAND "${Python3_EXECUTABLE}" -mpip install --no-deps -t "${NUMPY_PACKAGE}" numpy==${RFC_NUMPY_VERSION}
                OUTPUT_VARIABLE pip_output
                RESULT_VARIABLE pip_failed
        )
        if (pip_failed)
            message(FATAL_ERROR
                    "[Python] pip failed:\n"
                    "${pip_output}"
            )
        endif ()
    endif ()
    # Assemble new site-packages path (prefetch build-site numpy installation)
    set(tmp "${NUMPY_PACKAGE}")
    list(APPEND tmp "$ENV{PYTHONPATH}")
    list(JOIN tmp ${PATH_SEP} tmp)  # This works for linux and windows either
    set(ENV{PYTHONPATH} "${tmp}")
    message(STATUS
            "[Python] New modules path: "
            "$ENV{PYTHONPATH}"
    )
endif ()
# Need to use the NumPy C API
if (NOT Python3_NumPy_FOUND)
    message(STATUS "[Python] Didn't find Numpy, trying alternative search...")
    execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
            OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIRS
            RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )
    if (Python3_NumPy_NOTFOUND)
        message(FATAL_ERROR "[Python] NumPy headers not found")
    else ()
        message(STATUS "[Python] Numpy found.")
    endif ()
endif ()
# Get the active numpy version to be able to validate
execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.__version__, end='')"
        OUTPUT_VARIABLE Python3_NumPy_VERSION
        RESULT_VARIABLE Python3_NumPy_NOTFOUND
)
message(STATUS "[Python] Executable: ${Python3_EXECUTABLE}")
message(STATUS "[Python] Include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "[Python] Library: ${Python3_LIBRARIES}")
message(STATUS "[Python] Numpy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "[Python] Numpy version: ${Python3_NumPy_VERSION}")
# Build tag
if (NOT ("${RFC_NUMPY_VERSION}" STREQUAL ""))
    # Assemble build tags
    string(REPLACE . _ Numpy_tag ${Python3_NumPy_VERSION})
    set(Python_tag "cp${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}-${PYTHON_PLATFORM}")
    set(prebuild_name ${PROJECT_NAME}_npy_${Numpy_tag}.${Python_tag})
    message(STATUS "Set prebuild target name to \"${prebuild_name}\"")
else ()
    set(prebuild_name "")
endif ()

add_library(rfcnt SHARED
        src/rfcnt.cpp
        ${rfc_core_SOURCE_DIR}/rainflow.c
        ${rfc_core_SOURCE_DIR}/rainflow.h
        ${rfc_core_SOURCE_DIR}/config.h
)
target_include_directories(rfcnt PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
        ${rfc_core_SOURCE_DIR}
)
target_link_libraries(rfcnt PRIVATE ${Python3_LIBRARIES})
if (MINGW)
    # See: https://github.com/cython/cython/issues/3213
    target_link_options(rfcnt PRIVATE -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive)
    set_target_properties(rfcnt PROPERTIES PREFIX "")
endif ()
# Ignore rainflow configuration file and use fixed options:
target_compile_definitions(rfcnt PRIVATE
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
        # RFC_EXPORT_MEX
        RFC_UNIT_TEST
)
if (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set_target_properties(rfcnt PROPERTIES SUFFIX ".pyd")
endif ()
# Copy target
if ("${prebuild_name}" STREQUAL "")
    add_custom_command(TARGET rfcnt POST_BUILD COMMAND
            ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rfcnt> ${CMAKE_CURRENT_LIST_DIR}
            COMMENT "copying ...")
    install(TARGETS rfcnt LIBRARY DESTINATION ${CMAKE_CURRENT_LIST_DIR})
elseif (WIN32)
    set_target_properties(rfcnt PROPERTIES OUTPUT_NAME ${prebuild_name})
    add_custom_command(TARGET rfcnt POST_BUILD COMMAND
            ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rfcnt> ${CMAKE_CURRENT_LIST_DIR}/_ext
            COMMENT "copying to ./_ext ...")
endif ()
