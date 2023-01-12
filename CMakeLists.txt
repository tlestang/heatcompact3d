cmake_minimum_required(VERSION 3.10)
project(heatcompact3d LANGUAGES Fortran)
enable_testing()

file(GLOB SRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src/*.f90
)
add_library(heatcompact3d STATIC ${SRC})

file(GLOB TESTSRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/tests/*.f90
)
foreach(testsrc IN LISTS TESTSRC)
  get_filename_component(test_name ${testsrc} NAME_WE)
  add_executable(${test_name} ${testsrc})
  target_link_libraries(${test_name} PRIVATE heatcompact3d)

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
