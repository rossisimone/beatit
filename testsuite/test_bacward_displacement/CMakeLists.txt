SET(TESTNAMEV test_bacward_displacement)
message("=== Adding Test : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_bacward_displacement")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)


SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )

add_test(${TESTNAMEV} mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_bacward_displacement -i data.beat -pc_type lu)


SET(TESTNAMEV test_bacward_displacement_multimaterial)
message("=== Adding Test : ${TESTNAMEV}")

add_executable(${TESTNAMEV} multimaterial.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_bacward_displacement_multimaterial")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)



