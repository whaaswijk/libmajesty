cmake_minimum_required(VERSION 2.8)

# Embedded MiniSat directory
add_subdirectory(minisat)

set(SOURCE_FILES xmg.cpp mlputil.c  mlpext.cpp strashmap.cpp 
	mig_interface.cpp)