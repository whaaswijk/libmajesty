cmake_minimum_required(VERSION 2.8)

# Embedded MiniSat directory
add_subdirectory(minisat)

set(SOURCE_FILES xmg.cpp mlputil.c  mlpext.cpp strashmap.cpp 
	mig_interface.cpp convert.cpp bitset_utils.cpp truth_table_utils.cpp npn_canonization.cpp conversion_utils.cpp
	mig_functional_hashing_constants.cpp cut.cpp lut_cover.cpp lut_optimize.cpp maj_io.cpp)

add_library(libmajesty STATIC ${SOURCE_FILES})
target_link_libraries(libmajesty LINK_PUBLIC MiniSat ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
