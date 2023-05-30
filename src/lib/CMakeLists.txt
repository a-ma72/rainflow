# Rainflow core library

add_library(rfc_core STATIC ${CMAKE_CURRENT_LIST_DIR}/rainflow.c ${CMAKE_CURRENT_LIST_DIR}/rainflow.h ${CMAKE_CURRENT_LIST_DIR}/config.h)
target_link_libraries(rfc_core ${LIBM_LIBRARY})
target_include_directories(rfc_core PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_compile_definitions(rfc_core PRIVATE -DRFC_HAVE_CONFIG_H)

# Update sources for Python module
add_custom_command(TARGET rfc_core
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/rainflow.c ${CMAKE_CURRENT_LIST_DIR}/../python/src/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/rainflow.h ${CMAKE_CURRENT_LIST_DIR}/../python/src/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/rainflow.hpp ${CMAKE_CURRENT_LIST_DIR}/../python/src/
                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/config.h ${CMAKE_CURRENT_LIST_DIR}/../python/src/config.h_
                   COMMENT "Copy sources to python submodule" )
