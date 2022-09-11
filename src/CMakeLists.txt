set(ENSIGN_SOURCES ../src/lr/lr.cpp
                  ../src/lr/coefficients.cpp
                  ../src/generic/utility.cpp
                  ../src/generic/matrix.cpp
                  ../src/generic/timer.cpp
                  ../src/generic/fft.cpp
                  ../src/generic/netcdf.cpp)

set(ENSIGN_SOURCES_GPU ../src/lr/lr.cu
                  ../src/lr/coefficients.cu
                  ../src/generic/utility.cu
                  ../src/generic/matrix.cu
                  ../src/generic/timer.cu
                  ../src/generic/fft.cu
                  ../src/generic/kernels.cu
                  ../src/generic/netcdf.cpp)
