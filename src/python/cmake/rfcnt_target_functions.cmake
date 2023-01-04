# ==================================================== Functions ======================================================

# Function to get the NumPy C API version.
# The function takes one argument:
# - out_var: The variable name to store the detected NumPy C API version.

function(numpy_get_capi_version  out_var)
    # Initialize the output variable with "NOTFOUND".
    set(${out_var} "NOTFOUND")

    # Execute a Python command to get the NumPy include directory.
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include(), end='')"
        OUTPUT_VARIABLE Python3_NumPy_INCLUDE_DIR
        RESULT_VARIABLE Python3_NumPy_NOTFOUND
    )

    # Recursive search for the _numpyconfig.h file in the NumPy include directory.
    file(GLOB_RECURSE NUMPYCONFIG_H "${Python3_NumPy_INCLUDE_DIR}/_numpyconfig.h")

    # Check if the _numpyconfig.h file was found.
    if(NUMPYCONFIG_H)
        # Read the contents of the _numpyconfig.h file.
        file(STRINGS ${NUMPYCONFIG_H} NUMPYCONFIG)

        # Iterate over each line of the file contents.
        foreach(line ${NUMPYCONFIG})
            # Check if the line contains the NPY_API_VERSION definition.
            string(REGEX MATCH ".*#define[ \t]+NPY_API_VERSION[ \t]+" NPY_VERSION_FOUND ${line})

            # If the line contains the NPY_API_VERSION definition, extract the version number.
            if(NPY_VERSION_FOUND)
                string(REGEX REPLACE ".*#define[ \t]+NPY_API_VERSION[ \t]+0[xX]([0-9a-fA-F]+).*$" "\\1" NPY_VERSION "${line}")

                # If the version number was successfully extracted, convert it from hexadecimal to decimal.
                if (NPY_VERSION)
                    math(EXPR NPY_VERSION_DEC "0x${NPY_VERSION}" OUTPUT_FORMAT DECIMAL)

                    # Set the output variable with the converted version number and break the loop.
                    set(${out_var} ${NPY_VERSION_DEC} PARENT_SCOPE)
                    break()
                endif()
            endif()
        endforeach()
    endif ()
endfunction()


# Define a function to generate the appropriate pip argument for installing NumPy.
# The function takes two arguments:
# - out_var: The variable name to store the generated pip argument.
# - numpy_version: The specified version of NumPy.

function(numpy_pip_arg  out_var  numpy_version)
    # Check if the specified version is "oldest-supported-numpy".
    if (numpy_version STREQUAL "oldest-supported-numpy")
        # If the version is "oldest-supported-numpy", set the PIP_ARG to "oldest-supported-numpy".
        set(PIP_ARG "oldest-supported-numpy")
    else ()
        # Use a regular expression to check if the numpy_version is a simple version number without any qualifiers (e.g., "1.18.5").
        string(REGEX MATCH "^[0-9.]+$" NO_QUALIFIER "${numpy_version}")

        # If the numpy_version is a simple version number without qualifiers.
        if (NO_QUALIFIER)
            # Set the PIP_ARG to "numpy==<numpy_version>" to specify an exact version.
            set(PIP_ARG "numpy==${numpy_version}")
        else ()
            # Otherwise, set the PIP_ARG to "numpy<numpy_version>" to handle version specifiers (e.g., ">1.18").
            set(PIP_ARG "numpy${numpy_version}")
        endif ()
    endif ()
    # Set the output variable with the generated PIP_ARG and make it available in the parent scope.
    set(${out_var} ${PIP_ARG} PARENT_SCOPE)
endfunction()


# Define a function to generate a prebuild target name for the NumPy package.
# The function takes three arguments:
# - out_var: The variable name to store the generated prebuild target name.
# - numpy_version: The specified version of NumPy.
# - capi_version: The C API version of NumPy.

function(prebuild_target_name  out_var  numpy_version  capi_version)
    # Replace dots in the NumPy version with underscores to form a tag.
    string(REPLACE . _ Numpy_tag ${numpy_version})

    # Create a Python tag using the major and minor Python version and platform.
    set(Python_tag "cp${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}-${PYTHON_PLATFORM}")

    # Assemble the prebuild name using project name, NumPy tag, and Python tag.
    set(prebuild_name ${PROJECT_NAME}_npy_${Numpy_tag}_api_v${capi_version}.${Python_tag})

    # Set the output variable for the prebuild name, passing it to the parent scope.
    set(${out_var} ${prebuild_name} PARENT_SCOPE)
endfunction()



# ================================================== C module rfcnt ===================================================
function (rfcnt_library  numpy_version  dependencies  OUTPUT_NAME  OUTPUT_DIRECTORY  rfcnt_target)
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

        if (numpy_version)
            # If no version is given, the module is built in the current environment.
            # If a NumPy version string is specified, handle the installation and uninstallation
            # of that specific NumPy version.
            add_custom_target(
                    numpy_${rfcnt_target}_uninstall
                    COMMAND "${Python3_EXECUTABLE}" -mpip uninstall -y oldest-supported-numpy numpy
                    COMMENT "Remove previous NumPy installations"
                    DEPENDS ${dependencies}
            )
            # Create a custom target to uninstall any existing NumPy versions,
            # ensuring main module has been built already.

            add_custom_target(
                    numpy_${rfcnt_target}_install
                    COMMAND "${Python3_EXECUTABLE}" -mpip install numpy==${numpy_version}
                    COMMAND "${Python3_EXECUTABLE}" ${CMAKE_CURRENT_SOURCE_DIR}/cmake/numpy_get_include.py
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src
                    COMMENT "Installing numpy ${numpy_version}"
                    DEPENDS numpy_${rfcnt_target}_uninstall
            )
            # Create a custom target to install the specified NumPy version and get its include directory.

            add_dependencies(
                    ${rfcnt_target}
                    numpy_${rfcnt_target}_install
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
