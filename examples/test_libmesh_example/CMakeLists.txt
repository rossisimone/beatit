SET(TESTNAMEV test_libmesh_example)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_libmesh_example.x")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(Mesh "${CMAKE_CURRENT_BINARY_DIR}/mesh.xda")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/mesh.xda IS_NEWER_THAN ${Mesh} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/mesh.xda  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/mesh.xda  IS_NEWER_THAN ${Mesh} )