################################################################################
# catch2

find_package(Catch2 QUIET)
if(NOT catch2_FOUND)
  FetchContent_Declare(Catch2
    GIT_REPOSITORY "https://github.com/catchorg/Catch2"
    GIT_TAG "v3.7.1"
    GIT_SHALLOW ON
    GIT_PROGRESS ON
    OVERRIDE_FIND_PACKAGE
  )
  FetchContent_MakeAvailable(Catch2)
  find_package(Catch2 REQUIRED)
endif()


################################################################################

add_executable(test_index_functions test_index_functions.cpp)
add_executable(test_matrix_functions test_matrix_functions.cpp
    ../src/grid_class.cpp
    ../src/tree_class.cpp
)
add_executable(test_orthogonalization test_orthogonalization.cpp 
    ../src/grid_class.cpp
    ../src/matrix.cpp
    ../src/tree_class.cpp
)
add_executable(test_tree_h1 test_tree_h1.cpp
    ../src/grid_class.cpp
    ../src/matrix.cpp
    ../src/tree_class.cpp
)
add_executable(test_tree_h2 test_tree_h2.cpp
    ../src/grid_class.cpp
    ../src/matrix.cpp
    ../src/tree_class.cpp
)

add_test(NAME index_functions COMMAND test_index_functions)
add_test(NAME matrix_functions COMMAND test_matrix_functions)
add_test(NAME orthogonalization COMMAND test_orthogonalization)
add_test(NAME tree_h1 COMMAND test_tree_h1)
add_test(NAME tree_h2 COMMAND test_tree_h2)

target_link_libraries(test_index_functions
  PRIVATE
  Ensign::Ensign
  Catch2::Catch2WithMain
)
target_link_libraries(test_matrix_functions 
  PRIVATE
  Ensign::Ensign
  Catch2::Catch2WithMain
)
target_link_libraries(test_orthogonalization
  PRIVATE
  Ensign::Ensign
  Catch2::Catch2WithMain
)
target_link_libraries(test_tree_h1
  PRIVATE
  Ensign::Ensign
  Catch2::Catch2WithMain
)
target_link_libraries(test_tree_h2
  PRIVATE
  Ensign::Ensign
  Catch2::Catch2WithMain
)

target_include_directories(test_index_functions
  PRIVATE
  ${HDR_DIR}
)
target_include_directories(test_matrix_functions
  PRIVATE
  ${HDR_DIR}
  ${NETCDF_INCLUDE_DIR}
)
target_include_directories(test_orthogonalization
  PRIVATE
  ${HDR_DIR}
  ${NETCDF_INCLUDE_DIR}
)
target_include_directories(test_tree_h1
  PRIVATE
  ${HDR_DIR}
  ${NETCDF_INCLUDE_DIR}
)
target_include_directories(test_tree_h2
  PRIVATE
  ${HDR_DIR}
  ${NETCDF_INCLUDE_DIR}
)
