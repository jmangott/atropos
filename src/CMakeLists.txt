#####################################
# Set up the kinetic-cme executable #
#####################################

add_executable(${EXENAME1} main_kinetic.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(${EXENAME1} Ensign)
set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

target_include_directories(${EXENAME1} PRIVATE ${HEADER_FILES})

target_link_libraries(${EXENAME1} ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES})

if(${MKL_ENABLED})
    target_link_libraries(${EXENAME1} ${MKL_LIBRARIES})
else()
    target_link_libraries(${EXENAME1} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()

##########################################
# Set up the hierarchical-cme executable #
##########################################

add_executable(${EXENAME2} main_hierarchical.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(${EXENAME2} Ensign)
set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

target_include_directories(${EXENAME2} PRIVATE ${HEADER_FILES})

target_link_libraries(${EXENAME2} ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES} Eigen3::Eigen)

if(${MKL_ENABLED})
    target_link_libraries(${EXENAME2} ${MKL_LIBRARIES})
else()
    target_link_libraries(${EXENAME2} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()