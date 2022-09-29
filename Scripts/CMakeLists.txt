########################################################
# CMake file for building SSbarProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")
########################################################

ADD_SUBDIRECTORY( src )
ADD_EXECUTABLE( main_scripts.exe main_scripts.cc )
TARGET_LINK_LIBRARIES( main_scripts.exe Scripts )
