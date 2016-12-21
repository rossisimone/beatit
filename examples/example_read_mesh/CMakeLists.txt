add_executable(test_mesh main.cpp)
message("=== Adding Example : test_mesh")

set_target_properties(test_mesh PROPERTIES  OUTPUT "test_mesh.x")


target_link_libraries(test_mesh beatit)
target_link_libraries(test_mesh ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(test_mesh ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(test_mesh PROPERTIES LINKER_LANGUAGE CXX)
