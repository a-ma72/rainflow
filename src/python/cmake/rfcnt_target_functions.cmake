# ==================================================== Functions ======================================================

if (POLICY CMP0177)
  cmake_policy(SET CMP0177 OLD)
endif()


# ================================================== C module rfcnt ===================================================
function (rfcnt_library  OUTPUT_NAME  OUTPUT_DIRECTORY  rfcnt_target)
    # Create a new target that builds a rfcnt module.

    if (NOT TARGET ${rfcnt_target})
        # Check if the target does not already exist.

        set(
                INCLUDE_DIRS
                ${Python3_INCLUDE_DIRS}
                ${rfc_core_SOURCE_DIR}
                ${Python3_NumPy_INCLUDE_DIRS}
                ${PYBIND11_INCLUDE_DIR}
        )
        # Set the include directories for the target, including Python3 include directories,
        # source directory for rfc_core, NumPy include directories, and pybind11 headers
        # (rfcnt.cpp is a pybind11 extension module).

        set(
                RFCNT_SOURCES
                src/rfcnt.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/lib/rainflow.c
                ${CMAKE_CURRENT_SOURCE_DIR}/lib/rainflow.h
                ${CMAKE_CURRENT_SOURCE_DIR}/lib/rainflow.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/lib/config.h
        )
        # Set the source files for the target, including rfcnt.cpp and various rainflow library files.

        add_library(${rfcnt_target} SHARED ${RFCNT_SOURCES})
        # Create a shared library target with the specified source files.

        set_target_properties(
                ${rfcnt_target}
                PROPERTIES
                CXX_STANDARD 17
                CXX_STANDARD_REQUIRED ON
                CXX_EXTENSIONS OFF
        )
        # rfcnt.cpp is a pybind11 module and needs std::optional / structured bindings
        # (C++17); this only affects the .cpp source, rainflow.c keeps its own default
        # C standard.

        target_link_libraries(${rfcnt_target} PRIVATE ${Python3_LIBRARIES})
        # Link the Python3 libraries to the target.

        set_target_properties(
                ${rfcnt_target}
                PROPERTIES
                PREFIX ""
        )
        # Python imports extension modules by their exact file name (no "lib" prefix).
        # CMake's default SHARED-library prefix is already "" on MSVC, but on Unix-like
        # platforms (including MinGW, which follows Unix conventions) it defaults to
        # "lib", which would produce "librfcnt.so"/"librfcnt.dylib" instead of
        # "rfcnt.so"/"rfcnt.dylib" -- so this must be cleared unconditionally, not just
        # for MinGW.

        if (MINGW)
            # If using MinGW, set specific link options to address compatibility issues.
            # See: https://github.com/cython/cython/issues/3213
            target_link_options(
                    ${rfcnt_target}
                    PRIVATE
                    -static-libgcc
                    -static-libstdc++
                    -Wl,-Bstatic,--whole-archive
                    -lwinpthread
            )
        endif ()

        # Ignore rainflow configuration file and use fixed options:
        target_compile_definitions(
                ${rfcnt_target}
                PRIVATE
                NPY_TARGET_VERSION=NPY_1_19_API_VERSION
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
                RFC_EXPORT_MEX=0
                RFC_EXPORT_PY
                RFC_UNIT_TEST
                NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
                NPY_TARGET_VERSION=NPY_1_19_API_VERSION
        )
        # Define preprocessor macros for the target, configuring various options for the rainflow counting library.
        # NPY_NO_DEPRECATED_API/NPY_TARGET_VERSION matter for rfcnt.cpp, which pulls in NumPy's
        # C API transitively via <pybind11/numpy.h>; harmless no-ops for the plain-C rainflow.c.

        if (WIN32)
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    SUFFIX ".pyd"
            )
        endif ()
        # Windows python modules have suffix ".pyd"

        if (OUTPUT_DIRECTORY)
            # Set per-config output dirs so multi-config generators (VS, Xcode) don't
            # append a config subfolder (Debug/, Release/, …). Without these, the .pyd
            # lands in e.g. src/python/Debug/ and Python can't find it.
            set(_out_dir "${CMAKE_CURRENT_SOURCE_DIR}/${OUTPUT_DIRECTORY}")
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    RUNTIME_OUTPUT_DIRECTORY                 "${_out_dir}"
                    RUNTIME_OUTPUT_DIRECTORY_DEBUG           "${_out_dir}"
                    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO  "${_out_dir}"
                    RUNTIME_OUTPUT_DIRECTORY_RELEASE         "${_out_dir}"
                    RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL      "${_out_dir}"
                    LIBRARY_OUTPUT_DIRECTORY                 "${_out_dir}"
                    LIBRARY_OUTPUT_DIRECTORY_DEBUG           "${_out_dir}"
                    LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO  "${_out_dir}"
                    LIBRARY_OUTPUT_DIRECTORY_RELEASE         "${_out_dir}"
                    LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL      "${_out_dir}"
            )
            install(
                    TARGETS ${rfcnt_target}
                    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/${OUTPUT_DIRECTORY}
            )
            # If an output directory is specified, set the runtime output directory and install
            # the target library to that directory.
        endif ()

        if (OUTPUT_NAME)
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    OUTPUT_NAME ${OUTPUT_NAME}  # New executable name (my_app.exe on Windows)
            )
            # If an output name is specified, set the output name of the target.
        endif ()

        target_include_directories(
                ${rfcnt_target}
                BEFORE
                PRIVATE
                ${INCLUDE_DIRS}
        )
        # Set the include directories for the target, ensuring they are included before any other directories.

    endif ()
endfunction()
