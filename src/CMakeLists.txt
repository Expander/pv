add_library(pv b0.cpp)
target_include_directories(pv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(run run.cpp)
target_link_libraries(run PUBLIC pv)
