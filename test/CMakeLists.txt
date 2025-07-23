add_executable(rfc_test rfc_test.c rfc_wrapper_simple.cpp rfc_wrapper_advanced.cpp)
target_link_libraries(rfc_test PRIVATE rfc_core greatest)
target_compile_definitions(rfc_test PRIVATE RFC_HAVE_CONFIG_H)

if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT rfc_test)
endif ()