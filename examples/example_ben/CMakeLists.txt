SET(TESTNAMEV example_ben)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_ben")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/nh_TRI3_stabilized_incompressible.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/nh_TRI3_stabilized_incompressible.beat  IS_NEWER_THAN ${GetPotFile} )





