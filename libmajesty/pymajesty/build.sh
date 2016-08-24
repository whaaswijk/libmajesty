#!/usr/bin/env bash

# Compile libmajesty
cd ..
#rm -r build
mkdir -p build
cd build
cmake ..
make
cd ..

# Compile pymig
rm pymajesty.cpython-*
cd pymajesty
python setup.py build_ext --build-lib ..
cd ..
