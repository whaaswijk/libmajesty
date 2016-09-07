cmake_minimum_required(VERSION 2.8)

# Embedded MiniSat directory
add_subdirectory(minisat)

set(CMAKE_LIBRARY_APTH ${CMAKE_LIBRARY_PATH} $ENV{HIREDIS_HOME})

set(SOURCE_FILES xmg.cpp mlputil.c mlpmemory.c  mlpext.cpp strashmap.cpp 
	mig_interface.cpp convert.cpp bitset_utils.cpp truth_table_utils.cpp npn_canonization.cpp conversion_utils.cpp
	mig_functional_hashing_constants.cpp cut.cpp lut_cover.cpp lut_optimize.cpp maj_io.cpp function_store.cpp)

find_package(Boost COMPONENTS filesystem system REQUIRED)
find_package(BISON)
find_package(FLEX)
BISON_TARGET(Parser parser.y ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp)
FLEX_TARGET(Lexer lexer.l ${CMAKE_CURRENT_BINARY_DIR}/lexer.cpp)
ADD_FLEX_BISON_DEPENDENCY(Lexer Parser)


add_library(libmajesty STATIC ${SOURCE_FILES} ${BISON_Parser_OUTPUTS} ${FLEX_Lexer_OUTPUTS})
target_link_libraries(libmajesty LINK_PUBLIC MiniSat hiredis ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
