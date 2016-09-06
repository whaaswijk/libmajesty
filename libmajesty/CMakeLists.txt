cmake_minimum_required(VERSION 2.8)

project(libmajesty)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -Wall -Wextra -std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

#find_package(BISON)
#find_package(FLEX)
find_package(Threads)

#BISON_TARGET(Parser parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp)
#FLEX_TARGET(Lexer lexer.l ${CMAKE_CURRENT_BINARY_DIR}/lexer.cpp)
#ADD_FLEX_BISON_DEPENDENCY(Lexer Parser)

# Add the ./src and ./include dirs to the include path 
include_directories("$ENV{HIREDIS_HOME}")
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("$ENV{LIBABC_HOME}/src")

# Add the binary dir to the include path so we can find the config header
include_directories("${PROJECT_BINARY_DIR}")

# Add the src directory
add_subdirectory(src)

