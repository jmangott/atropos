#################################################################
# Define a short hand for the source files and the Ensign files #
#################################################################

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

set(SRC_FILES ${SRC}/main.cpp
              ${SRC}/reaction_class.cpp
)

# Set up the executable
add_executable(kinetic-cme ${SRC_FILES} ${ENSIGN_FILES})

######################################
# Include header files and libraries #
######################################

# Include project-specific header files
target_include_directories(kinetic-cme PRIVATE ${INC})

# Include Ensing header files
target_include_directories(kinetic-cme PRIVATE ${ENSIGN_DIR}/include)

# Include BLAS and LAPACK
option("MKL_ENABLED" OFF)
if(NOT ${MKL_ENABLED})
  set(OPENBLAS_LIBRARIES "${ENSIGN_BUILD}/extern/OpenBLAS/opt/OpenBLAS/lib/libopenblas.a")
  target_include_directories(kinetic-cme PRIVATE "${ENSIGN_BUILD}/extern/OpenBLAS/opt/OpenBLAS/include")
endif()
find_package(Threads REQUIRED)

# Include fftw
target_include_directories(kinetic-cme PRIVATE "${ENSIGN_BUILD}/extern/fftw/include")
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

# Link the libraries
target_link_libraries(kinetic-cme ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")