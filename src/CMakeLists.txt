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
)

add_executable(kinetic-cme main.cpp ${ENSIGN_FILES})

target_link_libraries(kinetic-cme ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
                  