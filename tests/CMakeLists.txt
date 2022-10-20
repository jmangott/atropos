enable_testing()

###################
# Test for K-Step #
###################

add_executable(k_step k_step.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(k_step Ensign)
set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

add_test(NAME k_step COMMAND k_step)

target_include_directories(k_step PRIVATE "${ENSIGN_BUILD}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})

target_link_libraries(k_step ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")