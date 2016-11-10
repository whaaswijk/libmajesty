cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -Wall -Wextra -std=c++11 -DLIN64")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "../../build/src")

include_directories("$ENV{LIBMAJESTY_HOME}/include")
include_directories("$ENV{HIREDIS_HOME}")
include_directories("$ENV{LIBABC_HOME}/src")

find_package(Boost COMPONENTS filesystem system REQUIRED)
find_package(Threads)
link_directories(${LINK_DIRECTORIES} $ENV{LIBMAJESTY_HOME}/build/src)
link_directories(${LINK_DIRECTORIES} $ENV{LIBMAJESTY_HOME}/build/src/minisat)
link_directories(${LINK_DIRECTORIES} $ENV{LIBABC_HOME})
link_directories(${LINK_DIRECTORIES} $ENV{HIREDIS_HOME})

set(SOURCE_FILES main.cpp)

add_executable(majesty_test ${SOURCE_FILES})
target_link_libraries(majesty_test LINK_PUBLIC libmajesty.a MiniSat libabc.a readline dl libhiredis.a cryptominisat5.so ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
