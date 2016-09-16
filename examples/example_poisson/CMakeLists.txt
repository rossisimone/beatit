SET(TESTNAME example_poisson)
add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "example_poisson")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/poisson.pot")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/poisson.pot  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/poisson.pot  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/poisson.pot  IS_NEWER_THAN ${GetPotFile} )


SET(TESTNAME example_restart)
add_executable(${TESTNAME} main_restart.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "example_restart")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)


