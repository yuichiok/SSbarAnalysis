########################################################
# cmake file for building NTupleProcessor
# @author Yuichi Okugawa, IJCLab
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
########################################################

# project name
PROJECT( SSbarProcessor )
SET(CMAKE_CXX_STANDARD 17)

ADD_SUBDIRECTORY( NtupleProcessor )
ADD_SUBDIRECTORY( SSbarLibrary )