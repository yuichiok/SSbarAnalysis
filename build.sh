#!/bin/bash

if [[ -d build ]]
then
    rm -rf build/
    echo "rm -rf build/"
fi

mkdir build
cd build
cmake -D CMAKE_C_COMPILER=/usr/bin/gcc -D CMAKE_CXX_COMPILER=/usr/bin/g++ ..
make
cd -