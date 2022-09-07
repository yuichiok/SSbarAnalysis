#!/bin/bash

if [[ -d build ]]
then
    rm -rf build/
    echo "rm -rf build/"
fi

mkdir build
cd build
cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ ..
make VERBOSE=1
cd -