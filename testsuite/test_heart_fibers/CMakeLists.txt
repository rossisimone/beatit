SET(TESTNAME test_heart_fibers)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_heart_fibers")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/fibers.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/fibers.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/fibers.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/fibers.beat  IS_NEWER_THAN ${GetPotFile} )

add_test(${TESTNAME} mpirun -n 4 ${CMAKE_CURRENT_BINARY_DIR}/test_heart_fibers -i fibers.beat)
