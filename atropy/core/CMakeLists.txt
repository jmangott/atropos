cmake_minimum_required(VERSION 3.27 FATAL_ERROR)
project(atropy
        VERSION 0.0.1
        DESCRIPTION "Dynamical low-rank solver for the kinetic CME"
        HOMEPAGE_URL https://git.uibk.ac.at/c7021158/kinetic-cme
        LANGUAGES CXX)
enable_testing()

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules")
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile" CACHE STRING "" FORCE)


################################################################################
# General options
################################################################################

option("MKL_ENABLED" "Enable MKL for Ensign" OFF)
option("OPENMP" "Enable OpenMP" OFF)


################################################################################
# Compiler options
################################################################################

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Warn all")
    set(CMAKE_CXX_FLAGS_DEBUG   "-g -check all -traceback")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip -xHOST")

else(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -ggdb3 -Wall -D_GLIBCXX_DEBUG -fno-omit-frame-pointer -ftrapv")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native")
endif()


################################################################################
# External libraries
################################################################################

################################################################################
# Eigen

find_package(Eigen3 3.4 REQUIRED NO_MODULE)


################################################################################
# NetCDF

find_package(NetCDF REQUIRED)


################################################################################
# cxxopts

find_package(cxxopts QUIET)
if(NOT cxxopts_FOUND)
  FetchContent_Declare(cxxopts
    GIT_REPOSITORY "https://github.com/jarro2783/cxxopts.git"
    GIT_TAG "v3.2.1"
    GIT_SHALLOW ON
    GIT_PROGRESS ON
    OVERRIDE_FIND_PACKAGE
  )
  FetchContent_MakeAvailable(cxxopts)
endif()


################################################################################
# Ensign

find_package(Ensign QUIET)
if(NOT ensign_FOUND)
  message(STATUS "Did not find Ensign, so we download it instead.")
  FetchContent_Declare(Ensign
    GIT_REPOSITORY https://github.com/leinkemmer/Ensign.git
    GIT_TAG "development"
    GIT_SHALLOW ON
    GIT_PROGRESS ON
    OVERRIDE_FIND_PACKAGE
  )
  FetchContent_MakeAvailable(Ensign)
  find_package(Ensign REQUIRED)
endif()


################################################################################
# Create short-hands
################################################################################

# Create `input` and `output` directories
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/input)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

# Define short-hands for the src directory and the files
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SRC_FILES ${SRC}/bug_integrator.cpp
              ${SRC}/grid_class.cpp
              ${SRC}/matrix.cpp
              ${SRC}/print_functions.cpp
              ${SRC}/ps_integrator.cpp
              ${SRC}/tree_class.cpp)
set(HDR_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

################################################################################
# Set up the program
################################################################################

add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${SRC_FILES})

# Add SOURCE_ROOT as a preprocessor directive
target_compile_definitions(${PROJECT_NAME} PRIVATE "SOURCE_ROOT=${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(${PROJECT_NAME} PRIVATE ${HDR_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE Ensign)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_INCLUDE_DIR})

if(cxxopts_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE cxxopts::cxxopts)
else()
  target_link_libraries(${PROJECT_NAME} PRIVATE cxxopts)
endif()

if(BUILD_TESTING)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()
