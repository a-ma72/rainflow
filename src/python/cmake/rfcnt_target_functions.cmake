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
        )
        # Set the include directories for the target, including Python3 include directories,
        # source directory for rfc_core, and NumPy include directories.

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

        target_link_libraries(${rfcnt_target} PRIVATE ${Python3_LIBRARIES})
        # Link the Python3 libraries to the target.

        if (MINGW)
            # If using MinGW, set specific link options and properties to address compatibility issues.
            # See: https://github.com/cython/cython/issues/3213
            target_link_options(
                    ${rfcnt_target}
                    PRIVATE
                    -static-libgcc
                    -static-libstdc++
                    -Wl,-Bstatic,--whole-archive
                    -lwinpthread
            )
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    PREFIX ""
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
        )
        # Define preprocessor macros for the target, configuring various options for the rainflow counting library.

        if (WIN32)
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    SUFFIX ".pyd"
            )
        endif ()
        # Windows python modules have suffix ".pyd"

        if (OUTPUT_DIRECTORY)
            set_target_properties(
                    ${rfcnt_target}
                    PROPERTIES
                    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${OUTPUT_DIRECTORY}  # Output directory
                    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${OUTPUT_DIRECTORY}  # Output directory
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
