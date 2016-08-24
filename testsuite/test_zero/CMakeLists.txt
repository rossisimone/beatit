add_executable(test_zero main.cpp)

set_target_properties(test_zero PROPERTIES  OUTPUT "test_zero.x")


target_link_libraries(test_zero beatit)
target_link_libraries(test_zero ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(test_zero ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(test_zero PROPERTIES LINKER_LANGUAGE CXX)
