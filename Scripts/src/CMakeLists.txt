########################################################
# CMake file for building SSbarProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")
########################################################

FILE(GLOB SOURCE_FILES
     "*.cc"
)

INCLUDE_DIRECTORIES(Scripts
  PUBLIC
  ../include
  ../../SSbarLibrary/include
)


ROOT_GENERATE_DICTIONARY(RootLibDict2 ../../SSbarLibrary/include/TreeStructures.hh LINKDEF ../LinkDef.hh)

ADD_LIBRARY(Scripts
  SHARED
  ${SOURCE_FILES}
  RootLibDict2
)

TARGET_INCLUDE_DIRECTORIES(Scripts
  PUBLIC
  ../include
  ../../SSbarLibrary/include
)

target_link_libraries(Scripts ${ROOT_LIBRARIES})

# install(TARGETS Scripts LIBRARY DESTINATION ../lib)
install(FILES ../../SSbarLibrary/include/TreeStructures.hh DESTINATION ../include)


# add_library(NtupleProcessor SHARED src/RootClass1.cpp src/RootClass2.cpp RootLibDict)