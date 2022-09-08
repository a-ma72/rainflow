project(rfcnt)
cmake_minimum_required(VERSION 3.12)

add_compile_definitions(
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

# Find Python 3
set(PYFTC_PYTHON_Version 3.8)
find_package(PythonInterp ${PYFTC_PYTHON_Version} REQUIRED)
find_package(PythonLibs ${PYFTC_PYTHON_Version} REQUIRED)
if( NOT Python_NumPy_FOUND )
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
message( STATUS "Python: ${PYTHON_EXECUTABLE}" )
message( STATUS "Python include: ${PYTHON_INCLUDE_DIR}" )
message( STATUS "Python lib: ${PYTHON_LIBRARY}" )
message( STATUS "Numpy: ${Python_NumPy_INCLUDE_DIRS}" )


add_library(rfcnt SHARED ../src/rainflow.c src/rfcnt.cpp)
set_target_properties(rfcnt PROPERTIES SUFFIX ".pyd")
target_include_directories(rfcnt PRIVATE ../src ${PYTHON_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(rfcnt PRIVATE ${PYTHON_LIBRARIES})
