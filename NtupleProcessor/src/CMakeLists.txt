# ADD_SUBDIRECTORY( NtupleProcessor )
# ADD_SUBDIRECTORY( SSbarLibrary )

CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)
SET(CMAKE_OSX_ARCHITECTURES "x86_64")

### DEPENDENCIES ############################################################

FIND_PACKAGE( ROOT REQUIRED )
SET(CMAKE_CXX_STANDARD 17)


FOREACH( pkg ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-c -Wall -fdiagnostics-color -std=gnu++17" )

# INCLUDE_DIRECTORIES( ${SSbarProcessor_SOURCE_DIR}/SSbarLibrary )
# LINK_DIRECTORIES( ${SSbarProcessor_SOURCE_DIR}/SSbarLibrary )
ADD_EXECUTABLE( SSbarProcessor main.cc )
INCLUDE_DIRECTORIES( ${SSbarProcessor_SOURCE_DIR}/NtupleProcessor/include )
TARGET_LINK_LIBRARIES(SSbarProcessor PUBLIC SSbarLibrary)

TARGET_INCLUDE_DIRECTORIES(SSbarProcessor PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/SSbarLibrary"
                          )