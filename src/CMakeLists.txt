add_subdirectory(lib)

add_library(greatest INTERFACE)
target_include_directories(greatest INTERFACE ${CMAKE_CURRENT_LIST_DIR}/3rd_party/greatest)

if(${RFC_EXPORT_MEX})
	add_subdirectory(matlab)
endif()
if(${RFC_EXPORT_PY})
	add_subdirectory(python)
endif()
