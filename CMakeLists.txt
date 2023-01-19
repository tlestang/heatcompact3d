cmake_minimum_required(VERSION 3.10)
project(heatcompact3d LANGUAGES Fortran)
enable_testing()

file(GLOB SRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src/*.f90
)

file(GLOB TESTSRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/tests/*.f90
)

file(GLOB CUDASRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src/cuda/*.f90
)

file(GLOB CUDATESTSRC
  LIST_DIRECTORIES false
  RELATIVE ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/tests/cuda/*.f90
)

# If the Fortran compiler is nvfortran, then append CUDA sources to
# global list of sources
string(REGEX MATCH "nvfortran" CUDAFORTRAN ${CMAKE_Fortran_COMPILER})
message("IS ${CUDAFORTRAN}")
if(CUDAFORTRAN)
  list(APPEND SRC ${CUDASRC})
  list(APPEND TESTSRC ${CUDATESTSRC})
endif()

message(${SRC})
add_library(heatcompact3d STATIC ${SRC})

# Build the test executables and register them to CTest.
# TODO Extract test config code into subdirectory, see
# https://github.com/tlestang/heatcompact3d/issues/8
foreach(testsrc IN LISTS TESTSRC)
  get_filename_component(test_name ${testsrc} NAME_WE)
  add_executable(${test_name} ${testsrc})
  target_link_libraries(${test_name} PRIVATE heatcompact3d)

  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
