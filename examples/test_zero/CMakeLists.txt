SET(TESTNAMEV test_zero)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "test_zero.x")


target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(${TESTNAMEV} ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(${TESTNAMEV} PROPERTIES LINKER_LANGUAGE CXX)
