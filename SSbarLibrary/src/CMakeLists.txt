########################################################
# CMake file for building SSbarProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")
########################################################

FILE(GLOB SOURCE_FILES
     "*.cc"
)

INCLUDE_DIRECTORIES(SSbarLibrary
  PUBLIC
  ../include
  ../../SSbarLibrary/include
)

ADD_LIBRARY(SSbarLibrary
  SHARED
  ${SOURCE_FILES}
)
