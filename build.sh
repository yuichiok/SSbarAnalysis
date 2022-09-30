#!/bin/bash

if [[ -d "build" ]]
then
    rm -rf build/
    echo "rm -rf build/"
fi

if compgen -G "*.exe" > /dev/null
then
    rm -rf *.exe
    echo "rm -rf *.exe"
fi

mkdir build
cd build

cmake ..
# cmake --build . --target Scripts

# make VERBOSE=1
make
cd -