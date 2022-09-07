########################################################
# cmake file for building NTupleProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
########################################################

# INCLUDE_DIRECTORIES( ${SSbarProcessor_SOURCE_DIR}/SSbarLibrary )
ADD_LIBRARY(SSbarLibrary timestamp.hh)