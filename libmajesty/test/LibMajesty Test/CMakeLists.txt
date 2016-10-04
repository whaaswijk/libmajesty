cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -Wall -Wextra -std=c++11")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "../../build/src")

include_directories("$ENV{LIBMAJESTY_HOME}/include")

find_package(Boost COMPONENTS filesystem system REQUIRED)
link_directories(${LINK_DIRECTORIES} $ENV{LIBMAJESTY_HOME}/build/src)
link_directories(${LINK_DIRECTORIES} $ENV{LIBMAJESTY_HOME}/build/src/minisat)

set(SOURCE_FILES main.cpp)

add_executable(majesty_test ${SOURCE_FILES})
target_link_libraries(majesty_test LINK_PUBLIC libmajesty MiniSat ${Boost_LIBRARIES})
