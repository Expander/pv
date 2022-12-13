add_definitions(-DTEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

add_subdirectory(alt)

function(add_pv_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} alt pv)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_pv_test(bench_b)
add_pv_test(test_a)
add_pv_test(test_b)
