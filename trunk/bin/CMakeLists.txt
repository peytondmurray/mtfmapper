file(GLOB sources ../src/*${cpp_ext})

add_executable(mtfmapper_bin ${sources})
set_target_properties(mtfmapper_bin PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_mapper)
target_link_libraries(mtfmapper_bin ${OpenCV_LIBS} ${FFTW_LIBRARIES} ${TBB_LIBRARIES})

