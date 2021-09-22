enable_testing()
add_executable(test_b test_b.cpp)
target_link_libraries(test_b pv)
add_test(NAME test_b COMMAND test_b)
