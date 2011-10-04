file(GLOB sources ../src/*${cpp_ext})

add_executable(mtfmapper_bin ${sources})
set_target_properties(mtfmapper_bin PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_mapper)
target_link_libraries(mtfmapper_bin ${OpenCV_LIBS} ${FFTW_LIBRARIES} ${TBB_LIBRARIES})

file(GLOB gen_sources ../src/generator/*${cpp_ext})
add_executable(generate_rectangle ${gen_sources})
set_target_properties(generate_rectangle PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_generate_rectangle)
target_link_libraries(generate_rectangle ${OpenCV_LIBS} ${TBB_LIBRARIES})

file(GLOB chart_sources ../src/test_chart_generator/*${cpp_ext})
add_executable(generate_test_chart ${chart_sources})
set_target_properties(generate_test_chart PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_generate_test_chart)
target_link_libraries(generate_test_chart ${OpenCV_LIBS} )

IF (QT4_FOUND)
  include(${QT_USE_FILE})
  file(GLOB qt_SRCS ${PROJECT_SOURCE_DIR}/src/gui/*${cpp_ext})
  qt4_automoc(${qt_SRCS})
  add_executable(mtf_mapper_gui ${qt_SRCS})
  target_link_libraries(mtf_mapper_gui ${QT_LIBRARIES})
  install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/mtf_mapper_gui DESTINATION bin)
ENDIF (QT4_FOUND)

install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/mtf_mapper ${PROJECT_BINARY_DIR}/bin/mtf_generate_rectangle ${PROJECT_BINARY_DIR}/bin/mtf_generate_test_chart DESTINATION bin)
