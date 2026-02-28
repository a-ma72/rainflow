# Rainflow core library

add_library(rfc_core STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/rainflow.c
    ${CMAKE_CURRENT_SOURCE_DIR}/rainflow.h
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h
)
target_link_libraries(rfc_core ${LIBM_LIBRARY})
target_include_directories(rfc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(rfc_core PRIVATE RFC_HAVE_CONFIG_H)
