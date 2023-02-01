cmake_minimum_required(VERSION 3.15.1)
project(kinetic-cme
        VERSION 0.0.1
        DESCRIPTION "Dynamical low-rank solver for the kinetic CME"
        HOMEPAGE_URL https://git.uibk.ac.at/c7021158/kinetic-cme
        LANGUAGES CXX)
enable_testing()

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
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
    set(CMAKE_CXX_FLAGS_PROFILE "-g1 -O3 -fno-inline -fno-optimize-sibling-calls")
endif()

##############################################################
# Define the actual files and folders that make up the build #
##############################################################

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(INC ${CMAKE_SOURCE_DIR}/include)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

############################################
# Find external packages and set libraries #
############################################

# Add Ensign
ExternalProject_Add(Ensign
  PREFIX "${EXTERNALS_DIR}"
  GIT_REPOSITORY https://github.com/leinkemmer/Ensign.git
  GIT_TAG "development"
  INSTALL_COMMAND ""
  UPDATE_COMMAND "" # do not rebuild external projects as part of the normal workflow
)
set(ENSIGN_SRC_DIR ${EXTERNALS_DIR}/src/Ensign)
set(ENSIGN_BUILD_DIR ${EXTERNALS_DIR}/src/Ensign-build)

# Find BLAS and LAPACK
option("MKL_ENABLED" OFF)
if(NOT ${MKL_ENABLED})
  set(OPENBLAS_LIBRARIES "${ENSIGN_BUILD_DIR}/extern/OpenBLAS/opt/OpenBLAS/lib/libopenblas.a")
endif()
find_package(Threads REQUIRED)

# Set fftw
set(FFTW_LIBRARIES "${ENSIGN_BUILD_DIR}/extern/fftw/lib/libfftw3.a" "${ENSIGN_BUILD_DIR}/extern/fftw/lib/libfftw3_omp.a")

# Find NetCDF
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

set(ENSIGN_FILES ${ENSIGN_SRC_DIR}/src/lr/lr.cpp
                 ${ENSIGN_SRC_DIR}/src/lr/coefficients.cpp
                 ${ENSIGN_SRC_DIR}/src/generic/utility.cpp
                 ${ENSIGN_SRC_DIR}/src/generic/matrix.cpp
                 ${ENSIGN_SRC_DIR}/src/generic/timer.cpp
                 ${ENSIGN_SRC_DIR}/src/generic/fft.cpp
                 ${ENSIGN_SRC_DIR}/src/generic/netcdf.cpp)

set(ENSIGN_FILES_GPU ${ENSIGN_SRC_DIR}/src/lr/lr.cu
                     ${ENSIGN_SRC_DIR}/src/lr/coefficients.cu
                     ${ENSIGN_SRC_DIR}/src/generic/utility.cu
                     ${ENSIGN_SRC_DIR}/src/generic/matrix.cu
                     ${ENSIGN_SRC_DIR}/src/generic/timer.cu
                     ${ENSIGN_SRC_DIR}/src/generic/fft.cu
                     ${ENSIGN_SRC_DIR}/src/generic/kernels.cu
                     ${ENSIGN_SRC_DIR}/src/generic/netcdf.cpp)

set(SRC_FILES ${SRC}/grid_class.cpp
              ${SRC}/index_functions.cpp
              ${SRC}/io_functions.cpp
              ${SRC}/print_functions.cpp
              ${SRC}/reaction_class.cpp
              ${SRC}/s_step_functions.cpp
              ${SRC}/weight_functions.cpp)

########################
# Include header files #
########################

set(HEADER_FILES ${INC}
                 ${ENSIGN_SRC_DIR}/include
                 ${ENSIGN_BUILD_DIR}/extern/OpenBLAS/opt/OpenBLAS/include
                 ${ENSIGN_BUILD_DIR}/extern/fftw/include)

###########################################
# Add src, test and output subdirectories #
###########################################

add_subdirectory(${SRC} ${BIN})
add_subdirectory(${CMAKE_SOURCE_DIR}/tests)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)
