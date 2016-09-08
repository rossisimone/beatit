SET(TESTNAME test_NumericVector)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_NumericVector")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})
#link_directories(${LIBMESH_DIR}/lib)


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

add_test(${TESTNAME} mpirun -n 2 ${CMAKE_CURRENT_BINARY_DIR}/test_NumericVector)
