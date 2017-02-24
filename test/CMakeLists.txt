cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -Wall -Wextra -std=c++11 -DLIN64")

include_directories(${HIREDIS_INCLUDE_DIR} ${ABC_INCLUDE_DIR})

find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )
find_package(Threads)
link_directories(${LINK_DIRECTORIES} ${PROJECT_BINARY_DIR}/lib)
link_directories(${LINK_DIRECTORIES} ${PROJECT_BINARY_DIR}/lib/minisat)
link_directories(${LINK_DIRECTORIES} $ENV{ABC_LIB_DIR})
link_directories(${LINK_DIRECTORIES} $ENV{HIREDIS_LIB_DIR})

add_executable(majesty_test main.cpp)
add_dependencies(majesty_test majesty)
target_link_libraries(majesty_test LINK_PUBLIC libmajesty.a MiniSat libabc readline curses dl libhiredis ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
