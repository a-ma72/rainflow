# Unit test

# cd build
# cmake -G "Visual Studio 16 2019 Win64" -DRFC_UNIT_TEST=1 ..
# cmake --build test --target rfc_test --config Release
# ctest -C Release

set(rfc_core_sources $<TARGET_PROPERTY:rfc_core,SOURCES>)
set(rfc_core_include_dir $<TARGET_PROPERTY:rfc_core,INCLUDE_DIRECTORIES>)

add_executable(rfc_test rfc_test.c rfc_wrapper_simple.cpp rfc_wrapper_advanced.cpp)
target_link_libraries(rfc_test PRIVATE rfc_core greatest)
target_compile_definitions(rfc_test PRIVATE -DRFC_HAVE_CONFIG_H)

target_include_directories(rfc_test PRIVATE greatest)

if(MSVC)
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT rfc_test )
endif()