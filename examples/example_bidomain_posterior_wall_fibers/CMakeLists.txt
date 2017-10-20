
SET(TESTNAME example_posterior_wall_fibers)
message("=== Adding Example : ${TESTNAME}")

add_executable(${TESTNAME} main_fibers.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "example_posterior_wall_fibers")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})
#link_directories(${LIBMESH_DIR}/lib)


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )






