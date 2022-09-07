CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")

project(SSbarProcessor CXX)
ADD_SUBDIRECTORY( src )

add_executable(main.exe main.cc)
target_link_libraries(main.exe NtupleProcessor)