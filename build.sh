#!/bin/bash

if [[ -d "build" ]]
then
    rm -rf build/
    echo "rm -rf build/"
fi

if compgen -G "main.exe" > /dev/null
then
    rm -rf main.exe
    echo "rm -rf main.exe"
fi

mkdir build
cd build

cmake ..
# cmake --build . --target Scripts

# make VERBOSE=1
make
cd -