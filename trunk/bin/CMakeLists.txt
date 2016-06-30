file(GLOB sources ../src/*${cpp_ext})

add_executable(mtfmapper_bin ${sources})
set_target_properties(mtfmapper_bin PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_mapper)
target_link_libraries(mtfmapper_bin ${OpenCV_LIBS})

file(GLOB gen_sources ../src/generator/*${cpp_ext})
add_executable(generate_rectangle ${gen_sources})
set_target_properties(generate_rectangle PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_generate_rectangle)
target_link_libraries(generate_rectangle ${OpenCV_LIBS})

file(GLOB chart_sources ../src/test_chart_generator/*${cpp_ext})
add_executable(generate_test_chart ${chart_sources})
set_target_properties(generate_test_chart PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME mtf_generate_test_chart)
target_link_libraries(generate_test_chart ${OpenCV_LIBS} )

# Qt5 
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Widgets)

get_target_property(q5name Qt5::Widgets LOCATION)

IF (EXISTS ${q5name})
  file(GLOB qt_SRCS ${PROJECT_SOURCE_DIR}/src/gui/*${cpp_ext})
  add_executable(mtf_mapper_gui ${qt_SRCS})
  qt5_use_modules(mtf_mapper_gui Widgets)
  target_link_libraries(mtf_mapper_gui ${QT_LIBRARIES})
  install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/mtf_mapper_gui DESTINATION bin)
ELSE (EXISTS ${q5name})
  message("QT5Widgets not found, not configuring GUI")
ENDIF (EXISTS ${q5name})

install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/mtf_mapper ${PROJECT_BINARY_DIR}/bin/mtf_generate_rectangle ${PROJECT_BINARY_DIR}/bin/mtf_generate_test_chart DESTINATION bin)
