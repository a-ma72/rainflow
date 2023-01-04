# ==================================================== Functions ======================================================
function(numpy_package_name  out_version_string  out_package_name  out_prebuild_name  out_pip_specifier  numpy_version)
    # Function to get various details about the NumPy package.

    string(REGEX REPLACE "[=~]*([0-9.]+)" "\\1" VERSION_STRING ${numpy_version})
    # Extract the version number from the numpy_version string (removing any prefix like '=' or '~=')

    string(SUBSTRING ${numpy_version} 0 2 VERSION_QUALIFIER)
    # Get the first two characters of numpy_version to check for version qualifier

    if (VERSION_QUALIFIER STREQUAL "~=")
        # If the version qualifier is '~=', adjust the version string to exclude the patch version
        # (VERSION_STRING is 1.20 for example)
        string(REGEX REPLACE "(.*)\\.[0-9]+$" "\\1" VERSION_STRING2 ${VERSION_STRING})
    else ()
        set(VERSION_STRING2 ${VERSION_STRING})
        set(VERSION_QUALIFIER "==")
    endif ()

    set(NUMPY_PACKAGE ${CMAKE_CURRENT_SOURCE_DIR}/numpy/${VERSION_STRING2})
    # Set the NumPy package path using the version string

    set(${out_package_name} ${NUMPY_PACKAGE} PARENT_SCOPE)
    set(${out_version_string} ${VERSION_STRING} PARENT_SCOPE)
    # Set the output variables for the package name and version string, passing them to the parent scope

    if (out_prebuild_name)
        # If out_prebuild_name is specified, assemble the build tags
        string(REPLACE . _ Numpy_tag ${VERSION_STRING2})
        # Replace dots in the NumPy version with underscores to form a tag

        set(Python_tag "cp${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}-${PYTHON_PLATFORM}")
        # Create a Python tag using the major and minor Python version and platform

        set(prebuild_name ${PROJECT_NAME}_npy_${Numpy_tag}.${Python_tag})
        # Assemble the prebuild name using project name, NumPy tag, and Python tag

        set(${out_prebuild_name} ${prebuild_name} PARENT_SCOPE)
        # Set the output variable for the prebuild name, passing it to the parent scope
    endif ()

    if (out_pip_specifier)
        # If out_pip_qualifier is specified, assemble the version qualifier for installation
        set(${out_pip_specifier} "${VERSION_QUALIFIER}${VERSION_STRING}" PARENT_SCOPE)
    endif ()
endfunction()


# ================================================== C module rfcnt ===================================================
function (rfcnt_library  numpy_version_string  dependencies  OUTPUT_NAME  OUTPUT_DIRECTORY  rfcnt_target)
    # Create a new target that builds a rfcnt module (or prebuild) with a specified NumPy version installed
    # before build.

    if (NOT TARGET ${rfcnt_target})
        # Check if the target does not already exist.

        set(
                INCLUDE_DIRS
                ${Python3_INCLUDE_DIRS}
                ${rfc_core_SOURCE_DIR}
                ${NUMPY_INCLUDE_DIRECTORIES}
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

        if (numpy_version_string)
            # If no version is given, the module is built in the current environment.
            numpy_package_name(VERSION_STRING  ""  ""  PIP_NUMPY_SPECIFIER  ${numpy_version})
            # If a NumPy version string is specified, handle the installation and uninstallation
            # of that specific NumPy version.
            add_custom_target(
                    numpy_${VERSION_STRING}_uninstall
                    COMMAND "${Python3_EXECUTABLE}" -mpip uninstall -y oldest-supported-numpy numpy
                    COMMENT "Remove previous NumPy installations"
                    DEPENDS ${dependencies} rfcnt_module
            )
            # Create a custom target to uninstall any existing NumPy versions,
            # ensuring main module has been built already.

            add_custom_target(
                    numpy_${VERSION_STRING}_install
                    COMMAND "${Python3_EXECUTABLE}" -mpip install numpy${PIP_NUMPY_SPECIFIER}
                    COMMAND "${Python3_EXECUTABLE}" ${CMAKE_CURRENT_SOURCE_DIR}/cmake/numpy_get_include.py
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
                    COMMENT "Installing numpy${PIP_NUMPY_SPECIFIER}"
                    DEPENDS numpy_${VERSION_STRING}_uninstall
            )
            # Create a custom target to install the specified NumPy version and get its include directory.

            add_dependencies(
                    ${rfcnt_target}
                    numpy_${VERSION_STRING}_install
            )
            # Set the target to depend on the NumPy installation target.
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
