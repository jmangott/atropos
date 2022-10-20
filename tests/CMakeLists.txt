###################
# Test for K-Step #
###################

add_executable(test_k_step test_k_step.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(test_k_step Ensign)
add_executable(test_index_functions test_index_functions.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(test_index_functions Ensign)

set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

add_test(NAME k_step COMMAND test_k_step)
add_test(NAME index_functions COMMAND test_index_functions)

target_include_directories(test_k_step PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})
target_include_directories(test_index_functions PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})

target_link_libraries(test_k_step ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
target_link_libraries(test_index_functions ${FFTW_LIBRARIES} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")