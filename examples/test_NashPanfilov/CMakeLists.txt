SET(TESTNAMEV test_NashPanfilov)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_NashPanfilov.x")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
#link_directories(${LIBMESH_DIR}/lib)


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)


SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.pot")
SET(Mesh "${CMAKE_CURRENT_BINARY_DIR}/mesh1.e")

IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.pot  IS_NEWER_THAN ${GetPotFile} )
    
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/mesh1.e IS_NEWER_THAN ${Mesh} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/mesh1.e  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/mesh1.e  IS_NEWER_THAN ${Mesh} )
    
