SET(TESTNAME example_em)
message("=== Adding Example : ${TESTNAME}")

add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "example_em")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/em.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/em.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/em.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/em.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/electrophysiology.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/elasticity.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/activation.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  IS_NEWER_THAN ${GetPotFile} )
