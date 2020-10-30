project(libmajesty)
cmake_minimum_required(VERSION 2.8)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -Wall -Wextra -std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DLIN64")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -DLIN64")

externalproject_add(
	abc-project
	GIT_REPOSITORY "https://github.com/berkeley-abc/abc.git"
	PREFIX "abc"
	CONFIGURE_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_COMMAND pwd && make ABC_USE_PIC=1 -j8 libabc.a && make
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)

SET(ABC_INCLUDE_DIR "${PROJECT_BINARY_DIR}/abc/src/abc-project/src")
SET(ABC_LIB_DIR "${PROJECT_BINARY_DIR}/abc/src/abc-project") 

add_library(libabc SHARED IMPORTED)
set_target_properties(libabc PROPERTIES IMPORTED_LOCATION ${ABC_LIB_DIR}/libabc.a)

externalproject_add(
	hiredis-project
	GIT_REPOSITORY "https://github.com/redis/hiredis.git"
	PREFIX "hiredis"
	CONFIGURE_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_COMMAND make -j8 libhiredis.a
	BUILD_IN_SOURCE 1
	INSTALL_COMMAND ""
)

set(HIREDIS_INCLUDE_DIR "${PROJECT_BINARY_DIR}/hiredis/src/hiredis-project")
set(HIREDIS_LIB_DIR "${PROJECT_BINARY_DIR}/hiredis/src/hiredis-project")

add_library(libhiredis SHARED IMPORTED)
set_target_properties(libhiredis PROPERTIES IMPORTED_LOCATION ${HIREDIS_LIB_DIR}/libhiredis.a)

# Add the ./src and ./include dirs to the include path 
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/include")

# Add the binary dir to the include path so we can find the config header
include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(shell)
