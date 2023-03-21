#########################
# Set up the executable #
#########################

add_executable(kinetic-cme main.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(kinetic-cme Ensign)
set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

target_include_directories(kinetic-cme PRIVATE ${HEADER_FILES})

target_link_libraries(kinetic-cme ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES})

if(${MKL_ENABLED})
    target_link_libraries(kinetic-cme ${MKL_LIBRARIES})
else()
    target_link_libraries(kinetic-cme ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()