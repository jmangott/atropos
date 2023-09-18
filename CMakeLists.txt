cmake_minimum_required(VERSION 3.15.1)
project(kinetic-cme
        VERSION 0.0.1
        DESCRIPTION "Dynamical low-rank solver for the kinetic CME"
        HOMEPAGE_URL https://git.uibk.ac.at/c7021158/kinetic-cme
        LANGUAGES CXX)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
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

else(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -ggdb3 -Wall -D_GLIBCXX_DEBUG -fsanitize=address -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -ftrapv")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native")
endif()

##############################################################
# Define the actual files and folders that make up the build #
##############################################################

# Define some directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(INC ${CMAKE_SOURCE_DIR}/include)
set(BIN ${CMAKE_SOURCE_DIR}/bin)

##############################
# Define the exectuable name #
##############################

set(EXENAME1 kinetic-cme)
set(EXENAME2 hierarchical-cme)

############################################
# Find external packages and set libraries #
############################################

# Options
option("MKL_ENABLED" "Enable MKL" OFF)
option("OPENMP" "Enable OpenMP" OFF)

# Ensign
ExternalProject_Add(Ensign
  PREFIX "${EXTERNALS_DIR}"
  GIT_REPOSITORY https://github.com/leinkemmer/Ensign.git
  GIT_TAG "development"
  INSTALL_COMMAND ""
  CMAKE_ARGS -DMKL_ENABLED=${MKL_ENABLED} -DOPENMP=${OPENMP_ENABLED}
  UPDATE_COMMAND "" # do not rebuild external projects as part of the normal workflow
)
set(ENSIGN_SRC_DIR ${EXTERNALS_DIR}/src/Ensign)
set(ENSIGN_BUILD_DIR ${EXTERNALS_DIR}/src/Ensign-build)

# OpenMP
if(${OPENMP})
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set(OPENMP_DEFINES "-D__OPENMP__")
    add_definitions(${OPENMP_DEFINES})
  endif()
endif()

# MKL/OpenBLAS
if(${MKL_ENABLED})
  set(MKL_DEFINES    "-D__MKL__" "-DMKL_ILP64" "-m64")
  set(MKL_INCLUDEDIR "$ENV{MKLROOT}/include")
  set(MKL_LIBRARIES  "-L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
  add_definitions(${MKL_DEFINES})
  include_directories(${MKL_INCLUDEDIR})
else()
  set(OPENBLAS_LIBRARIES "${ENSIGN_BUILD_DIR}/extern/OpenBLAS/opt/OpenBLAS/lib/libopenblas.a")
endif()
find_package(Threads REQUIRED)

# fftw
set(FFTW_LIBRARIES "${ENSIGN_BUILD_DIR}/extern/fftw/lib/libfftw3.a" "${ENSIGN_BUILD_DIR}/extern/fftw/lib/libfftw3_omp.a")

# NetCDF
find_package(NetCDF)
if(${NETCDF_FOUND})
  include_directories(${NETCDF_INCLUDE_DIR})
  message("test: ${NETCDF_INCLUDE_DIR}")
  add_definitions("-D__NETCDF__")
else()
  message("NetCDF has _not_ been found. If you want to write snapshots to disk please set the path in CMake appropriately.")
  set(NETCDF_LIBRARIES "")
endif()

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
              ${SRC}/integrators.cpp
              ${SRC}/io_functions.cpp
              ${SRC}/print_functions.cpp
              ${SRC}/reaction_class.cpp
              ${SRC}/s_step_functions.cpp
              ${SRC}/timer_class.cpp
              ${SRC}/tree_class.cpp
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