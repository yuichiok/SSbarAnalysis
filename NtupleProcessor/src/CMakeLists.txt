# ADD_SUBDIRECTORY( NtupleProcessor )
# ADD_SUBDIRECTORY( SSbarLibrary )

CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")

ADD_LIBRARY(NtupleProcessor
  SHARED
    EventAnalyzer.cc
    NtupleProcessor.cc
    TreeIterator.cc
)

TARGET_INCLUDE_DIRECTORIES(NtupleProcessor
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)