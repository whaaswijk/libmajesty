cmake_minimum_required(VERSION 2.8)

find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )

link_directories(${LINK_DIRECTORIES} ${PROJECT_BINARY_DIR}/lib)
link_directories(${LINK_DIRECTORIES} ${PROJECT_BINARY_DIR}/lib/minisat)
link_directories(${LINK_DIRECTORIES} $ENV{ABC_LIB_DIR})
link_directories(${LINK_DIRECTORIES} $ENV{HIREDIS_LIB_DIR})

add_executable(majesty_shell majesty_shell.cpp shell_env.cpp commands.cpp)
target_link_libraries(majesty_shell LINK_PUBLIC libmajesty.a MiniSat libabc readline curses dl libhiredis ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
