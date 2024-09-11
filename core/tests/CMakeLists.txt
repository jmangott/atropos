####################
# Set up the tests #
####################

add_executable(test_index_functions test_index_functions.cpp)
add_executable(test_matrix_functions test_matrix_functions.cpp)
add_executable(test_orthogonalization test_orthogonalization.cpp)
add_executable(test_tree_h1 test_tree_h1.cpp)
add_executable(test_tree_h2 test_tree_h2.cpp)

add_test(NAME index_functions COMMAND test_index_functions)
add_test(NAME matrix_functions COMMAND test_matrix_functions)
add_test(NAME orthogonalization COMMAND test_orthogonalization)
add_test(NAME tree_h1 COMMAND test_tree_h1)
add_test(NAME tree_h2 COMMAND test_tree_h2)

target_include_directories(test_index_functions PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})
target_include_directories(test_matrix_functions PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})
target_include_directories(test_orthogonalization PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})
target_include_directories(test_tree_h1 PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})
target_include_directories(test_tree_h2 PRIVATE "${ENSIGN_BUILD_DIR}/extern/Catch2/src/Catch2/single_include" ${HEADER_FILES})

target_link_libraries(test_index_functions ${LIBNAME})
target_link_libraries(test_matrix_functions ${LIBNAME})
target_link_libraries(test_orthogonalization ${LIBNAME})
target_link_libraries(test_tree_h1 ${LIBNAME})
target_link_libraries(test_tree_h2 ${LIBNAME})

if(${MKL_ENABLED})
    target_link_libraries(test_index_functions ${MKL_LIBRARIES})
    target_link_libraries(test_matrix_functions ${MKL_LIBRARIES})
    target_link_libraries(test_orthogonalization ${MKL_LIBRARIES})
    target_link_libraries(test_tree_h1 ${MKL_LIBRARIES})
    target_link_libraries(test_tree_h2 ${MKL_LIBRARIES})
else()
    target_link_libraries(test_index_functions ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
    target_link_libraries(test_matrix_functions ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
    target_link_libraries(test_orthogonalization ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
    target_link_libraries(test_tree_h1 ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
    target_link_libraries(test_tree_h2 ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()