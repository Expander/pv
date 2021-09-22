enable_testing()

function(add_pv_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} pv)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

add_pv_test(test_b)
