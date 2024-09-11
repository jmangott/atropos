#########################
# Set up the executable #
#########################

add_executable(${EXENAME} main.cpp)

target_compile_definitions(${EXENAME} PRIVATE "SOURCE_ROOT=${CMAKE_SOURCE_DIR}")
target_include_directories(${EXENAME} PRIVATE ${HEADER_FILES})
target_include_directories(${EXENAME} PRIVATE "${ENSIGN_BUILD_DIR}/extern/cxxopts/src/cxxopts/include")
target_link_libraries(${EXENAME} ${LIBNAME} Eigen3::Eigen)

if(${MKL_ENABLED})
    target_link_libraries(${EXENAME} ${MKL_LIBRARIES})
else()
    target_link_libraries(${EXENAME} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()