add_library(alt b0_neg.cpp degrassi.c gpt.cpp softsusy.cpp)
target_include_directories(alt PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
