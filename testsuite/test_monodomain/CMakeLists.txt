SET(TESTNAME test_monodomain)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_monodomain")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nash_panfilov.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.pot  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nash_panfilov_amr.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.pot  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/ORd.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/ORd.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/ORd.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/ORd.pot  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/TP06.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/TP06.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/TP06.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/TP06.pot  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/Grandi11.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.pot  IS_NEWER_THAN ${GetPotFile} )



add_test(${TESTNAME}_NashPanfilov mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i nash_panfilov.pot)
add_test(${TESTNAME}_NashPanfilov_AMR mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i nash_panfilov_amr.pot)
add_test(${TESTNAME}_ORd mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i ORd.pot)
add_test(${TESTNAME}_TP06 mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i TP06.pot)
add_test(${TESTNAME}_Grandi11 mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i Grandi11.pot)