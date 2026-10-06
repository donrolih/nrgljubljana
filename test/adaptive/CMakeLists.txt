# Tests of adaptive broadening (Lee and Weichselbaum, PRB 94, 235127 (2016)). Each test has its own 'run' script and
# compares runs with each other or with closed forms; see README.md. QS is always compiled, so these always run.

find_package(Perl REQUIRED)

file(STRINGS simple_tests adaptive_tests)
foreach(test ${adaptive_tests})
  add_test(NAME ${test}
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${test}/run ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${test}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${TEST_BUILD_ENV}" LABELS adaptive)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
