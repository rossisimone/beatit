SET(TESTNAME test_fibers)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_fibers")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nash_panfilov.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  IS_NEWER_THAN ${GetPotFile} )

add_test(${TESTNAME} mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_fibers -i nash_panfilov.pot)
