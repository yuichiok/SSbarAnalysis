########################################################
# CMake file for building SSbarProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")
########################################################

FILE(GLOB SOURCE_FILES
     "*.cc"
)

INCLUDE_DIRECTORIES(NtupleProcessor
  PUBLIC
  ../include
  ../../SSbarLibrary/include
)


ROOT_GENERATE_DICTIONARY(RootLibDict ../../SSbarLibrary/include/TreeStructures.hh LINKDEF ../LinkDef.hh)

ADD_LIBRARY(NtupleProcessor
  SHARED
  ${SOURCE_FILES}
  RootLibDict
)

TARGET_INCLUDE_DIRECTORIES(NtupleProcessor
  PUBLIC
  ../include
  ../../SSbarLibrary/include
)

target_link_libraries(NtupleProcessor ${ROOT_LIBRARIES})

install(TARGETS NtupleProcessor LIBRARY DESTINATION ../lib)
install(FILES ../../SSbarLibrary/include/TreeStructures.hh DESTINATION ../include)


# add_library(NtupleProcessor SHARED src/RootClass1.cpp src/RootClass2.cpp RootLibDict)