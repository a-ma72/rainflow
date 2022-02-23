project(rfcnt)
cmake_minimum_required(VERSION 3.12)

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
message( STATUS "Numpy: ${Python_NumPy_INCLUDE_DIRS}" )


add_library(rfcnt SHARED ../src/rainflow.c src/rfcnt.cpp ${PYTHON_LIBRARY})
target_include_directories(rfcnt PRIVATE ../src ${PYTHON_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(rfcnt PRIVATE ${PYTHON_LIBRARIES})
