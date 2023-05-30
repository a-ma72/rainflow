project(rfcnt)
cmake_minimum_required(VERSION 3.14)


# Find Python 3
if(DEFINED ENV{VIRTUAL_ENV})
    set(Python_ROOT_DIR $ENV{VIRTUAL_ENV})
elseif(DEFINED ENV{PYTHONHOME})
    set(Python_ROOT_DIR $ENV{PYTHONHOME})
endif()
if(NOT ("${RFC_NUMPY_VERSION}" STREQUAL ""))
    exec_program(pip ARGS "install --no-deps numpy==${RFC_NUMPY_VERSION}" OUTPUT_VARIABLE pip_output)
    message(Test ${pip_output})
endif()
find_package(Python ${RFC_PYTHON_VERSION} EXACT REQUIRED)
find_package(PythonInterp ${RFC_PYTHON_VERSION} EXACT REQUIRED)
find_package(PythonLibs ${RFC_PYTHON_VERSION} EXACT REQUIRED)
if(NOT Python_NumPy_FOUND)
    # We also need to use the NumPy C API
    MESSAGE(STATUS ${PYTHON_EXECUTABLE})
    exec_program(${PYTHON_EXECUTABLE}
                 ARGS "-c \"import numpy; print(numpy.get_include())\""
                 OUTPUT_VARIABLE Python_NumPy_INCLUDE_DIRS
                 RETURN_VALUE Python_NumPy_NOTFOUND
                )
    if(Python_NumPy_NOTFOUND)
        message(FATAL_ERROR "NumPy headers not found")
    endif()
endif()
exec_program(${PYTHON_EXECUTABLE}
             ARGS "-c \"import numpy; print(numpy.__version__)\""
             OUTPUT_VARIABLE Python_NumPy_VERSION
             RETURN_VALUE Python_NumPy_NOTFOUND
            )
message( STATUS "Python: ${PYTHON_EXECUTABLE}" )
message( STATUS "Python include: ${PYTHON_INCLUDE_DIR}" )
message( STATUS "Python lib: ${PYTHON_LIBRARY}" )
message( STATUS "Numpy: ${Python_NumPy_INCLUDE_DIRS}" )
message( STATUS "Numpy version: ${Python_NumPy_VERSION}" )

if(NOT ("${RFC_NUMPY_VERSION}" STREQUAL "" ))
    string(REPLACE . _ Numpy_tag ${Python_NumPy_VERSION})
    set(Python_tag cp${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-win_amd64)
    set(prebuild_name ${PROJECT_NAME}_npy_${Numpy_tag}.${Python_tag})
else()
    set(prebuild_name "")
endif()

set(rfc_core_SOURCES $<TARGET_PROPERTY:rfc_core,SOURCES>)
set(rfc_core_INCLUDE_DIRS $<TARGET_PROPERTY:rfc_core,INCLUDE_DIRECTORIES>)
add_library(rfcnt SHARED src/rfcnt.cpp ${rfc_core_SOURCES})
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set_target_properties(rfcnt PROPERTIES SUFFIX ".pyd")
endif()
target_include_directories(rfcnt PRIVATE ../src ${PYTHON_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS} ${rfc_core_INCLUDE_DIRS})
target_link_libraries(rfcnt PRIVATE ${PYTHON_LIBRARIES})
if(MINGW)
    # See: https://github.com/cython/cython/issues/3213
    target_link_options(rfcnt PRIVATE -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive)
    set_target_properties(rfcnt PROPERTIES PREFIX "")
endif()
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
    # RFC_TEST
)
# Copy target
if("${prebuild_name}" STREQUAL "")
    add_custom_command(TARGET rfcnt POST_BUILD COMMAND 
        ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rfcnt> ${CMAKE_CURRENT_LIST_DIR}
        COMMENT "copying ...")
    install(TARGETS rfcnt LIBRARY DESTINATION ${CMAKE_CURRENT_LIST_DIR})
elseif(WIN32)
    set_target_properties(rfcnt PROPERTIES OUTPUT_NAME ${prebuild_name})
    add_custom_command(TARGET rfcnt POST_BUILD COMMAND 
        ${CMAKE_COMMAND} -E copy $<TARGET_FILE:rfcnt> ${CMAKE_CURRENT_LIST_DIR}/_ext
        COMMENT "copying to ./_ext ...")
endif()
