cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-logical-op-parentheses -std=c++11")

add_library(MiniSat STATIC Solver.cpp Options.cpp)

include_directories(".")

# Make sure the compiler can find include files for our library
# when other libraries or executables link to it
target_include_directories(MiniSat PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
