# Host tests for the standalone RF_* engine.
# test_rainflow.c supplies its own RF_DAMAGE_LUT — do not link rainflow_config.c.

add_executable(rf_embedded_test
    rainflow.c
    test/test_rainflow.c
)
target_include_directories(rf_embedded_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Self-contained reference vs O(1) flush; does not compile rainflow.c.
add_executable(rf_embedded_flush_test
    test/test_RF_FlushResiduumRepeated.c
)

add_test(NAME rf_embedded_test COMMAND rf_embedded_test)
add_test(NAME rf_embedded_flush_test COMMAND rf_embedded_flush_test)
