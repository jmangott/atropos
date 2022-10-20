#########################
# Set up the executable #
#########################

add_executable(kinetic-cme main.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(kinetic-cme Ensign)

target_include_directories(kinetic-cme PRIVATE ${HEADER_FILES})

target_link_libraries(kinetic-cme ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")