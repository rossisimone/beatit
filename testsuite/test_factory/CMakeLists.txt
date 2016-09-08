SET(TESTNAME test_factory)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_factory")

target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

add_test(${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/test_factory)