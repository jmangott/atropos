cmake_minimum_required(VERSION 3.15.1)
project(kinetic-cme
        VERSION 0.0.1
        DESCRIPTION "Dynamical low-rank solver for the kinetic CME"
        HOMEPAGE_URL https://git.uibk.ac.at/c7021158/kinetic-cme
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile" CACHE STRING "" FORCE)

include(ExternalProject)
set(EXTERNALS_DIR ${CMAKE_BINARY_DIR}/external)

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

##############
# Add Ensign #
##############

ExternalProject_Add(Ensign
  PREFIX "${EXTERNALS_DIR}"
  GIT_REPOSITORY https://github.com/leinkemmer/Ensign.git
  GIT_TAG "development"
  SOURCE_DIR "${EXTERNALS_DIR}/Ensign"
  BINARY_DIR "${EXTERNALS_DIR}/Ensign/build"
  CONFIGURE_COMMAND cmake ..
  INSTALL_COMMAND ""
  UPDATE_COMMAND "" # do not rebuild external projects as part of the normal workflow
)

set(ENSIGN_DIR ${EXTERNALS_DIR}/Ensign)
set(ENSIGN_BUILD ${EXTERNALS_DIR}/Ensign/build)
# add_subdirectory(${ENSIGN_DIR} ${ENSIGN_BUILD})

##############################################################
# Define the actual files and folders that make up the build #
##############################################################

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(INC ${CMAKE_SOURCE_DIR}/include)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

######################################
# Include header files and libraries #
######################################

# Include BLAS and LAPACK
option("MKL_ENABLED" OFF)
if(NOT ${MKL_ENABLED})
  set(OPENBLAS_LIBRARIES "${ENSIGN_BUILD}/extern/OpenBLAS/opt/OpenBLAS/lib/libopenblas.a")
endif()
find_package(Threads REQUIRED)

# Include fftw
set(FFTW_LIBRARIES "${ENSIGN_BUILD}/extern/fftw/lib/libfftw3.a" "${ENSIGN_BUILD}/extern/fftw/lib/libfftw3_omp.a")

# Include NetCDF
# find_package(NetCDF)
# if(${NETCDF_FOUND})
#   include_directories(${NETCDF_INCLUDE_DIR})
#   message("test: ${NETCDF_INCLUDE_DIR}")
#   add_definitions("-D__NETCDF__")
# else()
#   message("NetCDF has _not_ been found. If you want to write snapshots to disk please set the path in CMake appropriately.")
#   set(NETCDF_LIBRARIES "")
# endif()

###############################################################
# Define a short hand for source, the Ensign and header files #
###############################################################

set(ENSIGN_FILES ${ENSIGN_DIR}/src/lr/lr.cpp
                 ${ENSIGN_DIR}/src/lr/coefficients.cpp
                 ${ENSIGN_DIR}/src/generic/utility.cpp
                 ${ENSIGN_DIR}/src/generic/matrix.cpp
                 ${ENSIGN_DIR}/src/generic/timer.cpp
                 ${ENSIGN_DIR}/src/generic/fft.cpp
                 ${ENSIGN_DIR}/src/generic/netcdf.cpp)

set(ENSIGN_FILES_GPU ${ENSIGN_DIR}/src/lr/lr.cu
                     ${ENSIGN_DIR}/src/lr/coefficients.cu
                     ${ENSIGN_DIR}/src/generic/utility.cu
                     ${ENSIGN_DIR}/src/generic/matrix.cu
                     ${ENSIGN_DIR}/src/generic/timer.cu
                     ${ENSIGN_DIR}/src/generic/fft.cu
                     ${ENSIGN_DIR}/src/generic/kernels.cu
                     ${ENSIGN_DIR}/src/generic/netcdf.cpp)

set(SRC_FILES ${SRC}/reaction_class.cpp
              ${SRC}/index_functions.cpp
              ${SRC}/io_functions.cpp)

# add_library(libensign STATIC ${ENSIGN_FILES})

# ########################
# # Include header files #
# ########################

set(HEADER_FILES ${INC}
                 ${ENSIGN_DIR}/include
                 ${ENSIGN_BUILD}/extern/OpenBLAS/opt/OpenBLAS/include
                 ${ENSIGN_BUILD}/extern/fftw/include
)

# target_include_directories(libensign PRIVATE ${HEADER_FILES})

###################################
# Add src and test subdirectories #
###################################

add_subdirectory(${SRC} ${BIN})
add_subdirectory(${CMAKE_SOURCE_DIR}/tests)