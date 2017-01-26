#!/usr/bin/env bash

ORIGINAL_DIR="$( pwd )"
echo ${ORIGINAL_DIR}

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Compile libmajesty
cd ${SCRIPT_DIR}/..
#cd libmajesty/libmajesty
#rm -r build
mkdir -p build
cd build
cmake ..
make

# Compile pymig
cd ${SCRIPT_DIR}
rm -f ${ORIGINAL_DIR}/majespy.cpython-*
python setup.py build_ext --build-lib ${ORIGINAL_DIR}
cd ${ORIGINAL_DIR}
