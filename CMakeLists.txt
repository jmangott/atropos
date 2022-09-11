cmake_minimum_required(VERSION 3.15.1)
project(kinetic-cme
        VERSION 0.0.1
        DESCRIPTION "Dynamical low-rank solver for the kinetic CME"
        HOMEPAGE_URL https://git.uibk.ac.at/c7021158/kinetic-cme
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")

set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile" CACHE STRING "" FORCE)

####################
# Compiler options #
####################

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Warn all")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -check all -traceback")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip -xHOST")
    set(CMAKE_CXX_FLAGS_PROFILE "-debug minimal -O3 -fno-inline -no-ip -no-ipo -fno-omit-frame-pointer -fno-optimize-sibling-calls")

else(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
    set(CMAKE_CXX_FLAGS_PROFILE "-g1 -O3 -fno-inline -fno-optimize-sibling-calls")
endif()

################
# BLAS options #
################

option("MKL_ENABLED" OFF)
find_package(Threads REQUIRED)

##############
# Add Ensign #
##############

FetchContent_Declare(Ensign
  PREFIX "${PROJECT_BINARY_DIR}/extern/Ensign"
  GIT_REPOSITORY "https://github.com/leinkemmer/Ensign.git"
  GIT_TAG "development"
  # CONFIGURE_COMMAND cmake ../Ensign
  # BUILD_COMMAND make
  # INSTALL_COMMAND ""
  # UPDATE_COMMAND "" # do not rebuild external projects as part of the normal workflow
)
FetchContent_MakeAvailable(Ensign)
add_subdirectory(Ensign)
# set(ENSIGN_DIR ${PROJECT_BINARY_DIR}/extern/Ensign/src/Ensign/src)
# include_directories("${PROJECT_BINARY_DIR}/extern/Ensign/src/Ensign/include")

##############################################################
# Define the actual files and folders that make up the build #
##############################################################

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(INC ${CMAKE_SOURCE_DIR}/include)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

# Add directories
include_directories(${INC})
add_subdirectory(${SRC} ${BIN})
