SET(TESTNAME test_monodomain)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_monodomain")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nash_panfilov.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nash_panfilov_amr.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nash_panfilov_amr.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/ORd.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/ORd.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/ORd.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/ORd.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/TP06.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/TP06.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/TP06.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/TP06.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/Grandi11.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/Grandi11.beat  IS_NEWER_THAN ${GetPotFile} )



add_test(${TESTNAME}_NashPanfilov mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i nash_panfilov.beat)
add_test(${TESTNAME}_NashPanfilov_AMR mpirun -n 4 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i nash_panfilov_amr.beat)
add_test(${TESTNAME}_ORd mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i ORd.beat)
add_test(${TESTNAME}_TP06 mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i TP06.beat)
add_test(${TESTNAME}_Grandi11 mpirun -n 8 ${CMAKE_CURRENT_BINARY_DIR}/test_monodomain -i Grandi11.beat)