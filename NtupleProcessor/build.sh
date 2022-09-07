#!/bin/bash

if [[ -d build ]]
then
    rm -rf build/
    echo "rm -rf build/"
fi

mkdir build
cd build
# cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ ..
# cmake -DCMAKE_C_COMPILER=/pbs/software/centos-7-x86_64/gcc/7.3.0/bin/gcc -DCMAKE_CXX_COMPILER=/pbs/software/centos-7-x86_64/gcc/7.3.0/bin/g++ ..
cmake ..
make VERBOSE=1
cd -