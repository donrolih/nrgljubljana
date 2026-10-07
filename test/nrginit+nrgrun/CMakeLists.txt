# Regression tests

file(STRINGS simple_tests simple_tests)

foreach(test ${simple_tests})
  # A test with its own run script compares runs with each other and stores no reference.
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run)
    add_test(NAME ${test}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}" LABELS adaptive)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
    continue()
  endif()
  file(GLOB inputs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${test}/*.m)
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    add_chain_test(${test} ${backend} MODE initializer SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${test} INPUTS param ${inputs}
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtest ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test} ${Mathematica_KERNEL_EXECUTABLE})
  endforeach()
endforeach()
