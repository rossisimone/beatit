SET(TESTNAMEV test_quasistatic_elasticity)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_quasistatic_elasticity")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)


SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/linear_TRI3.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/linear_TRI3_stabilized.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/linear_TRI3_stabilized_incompressible.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized_incompressible.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/linear_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nh_TRI3.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nh_TRI3_stabilized.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nh_TRI3_stabilized_incompressible.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )

add_test(${TESTNAMEV}_linear_TRI3 mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i linear_TRI3.beat -pc_type lu)
add_test(${TESTNAMEV}_linear_TRI3_stabilized mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i linear_TRI3_stabilized.beat -pc_type lu)
add_test(${TESTNAMEV}_linear_TRI3_stabilized_incompressible mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i linear_TRI3_stabilized_incompressible.beat -pc_type lu)
add_test(${TESTNAMEV}_nh_TRI3 mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i nh_TRI3.beat -pc_type lu)
add_test(${TESTNAMEV}_nh_TRI3_stabilized mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i nh_TRI3_stabilized.beat -pc_type lu)
add_test(${TESTNAMEV}_nh_TRI3_stabilized_incompressible mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_quasistatic_elasticity -i nh_TRI3_stabilized_incompressible.beat -pc_type lu)





