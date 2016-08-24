#!/usr/bin/env bash

# Compile libmajesty
cd libmajesty/libmajesty
#rm -r build
mkdir -p build
cd build
cmake ..
make
cd ../../..

# Compile pymig
rm pymig.cpython-*
cd pymig
python setup.py build_ext --build-lib ..
cd ..
