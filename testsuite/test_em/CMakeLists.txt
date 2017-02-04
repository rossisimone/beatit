SET(TESTNAMEV test_em)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_em")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)


SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/activation.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/elasticity.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/electrophysiology.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/em.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/em.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/em.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/em.beat  IS_NEWER_THAN ${GetPotFile} )

add_test(${TESTNAMEV}_NashPanfilov mpirun -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_em -i em.beat -pc_type lu)





