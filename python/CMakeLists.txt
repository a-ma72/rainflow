project(pyrfc)

# Find Python 3
find_package(Python3 COMPONENTS Interpreter Development)

# We also need to use the NumPy C API
exec_program(${_Python3_EXECUTABLE}
             ARGS "-c \"import numpy; print(numpy.get_include())\""
             OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
             RETURN_VALUE NUMPY_NOT_FOUND
            )
if(NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif()

add_library(pyrfc SHARED pyrfc.cpp)
target_include_directories(pyrfc PRIVATE ${_Python3_INCLUDE_DIR})
target_link_libraries(pyrfc PRIVATE ${_Python3_LIBRARY_RELEASE})
