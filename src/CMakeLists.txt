#########################
# Set up the executable #
#########################

add_executable(${EXENAME} main.cpp ${SRC_FILES} ${ENSIGN_FILES})
add_dependencies(${EXENAME} Ensign)
set_source_files_properties(${ENSIGN_FILES} PROPERTIES GENERATED TRUE)

target_include_directories(${EXENAME} PRIVATE ${HEADER_FILES})

target_link_libraries(${EXENAME} ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES})

if(${MKL_ENABLED})
    target_link_libraries(${EXENAME} ${MKL_LIBRARIES})
else()
    target_link_libraries(${EXENAME} ${OPENBLAS_LIBRARIES} Threads::Threads "-lgfortran")
endif()