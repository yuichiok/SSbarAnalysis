########################################################
# CMake file for building SSbarProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")
########################################################

FILE(GLOB SOURCE_FILES
     "*.cc"
)

ADD_LIBRARY(NtupleProcessor
  SHARED
  ${SOURCE_FILES}
)

TARGET_INCLUDE_DIRECTORIES(NtupleProcessor
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/../SSbarLibrary
)