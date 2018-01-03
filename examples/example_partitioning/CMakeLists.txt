SET(TESTNAMEV example_partitioning)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_partitioning")


#target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(${TESTNAMEV} ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV example_restart_partitioned_mesh)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main_restart.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_restart_partitioned_mesh")


#target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(${TESTNAMEV} ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)
