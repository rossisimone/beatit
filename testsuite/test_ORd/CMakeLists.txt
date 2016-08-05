add_executable(test_ORd main.cpp)

set_target_properties(test_ORd PROPERTIES  OUTPUT "test_ORd.x")


target_link_libraries(test_ORd beatit)
include_directories ("${PROJECT_SOURCE_DIR}/src")
#CMAKE_SOURCE_DIR
#target_include_directories(test_zero ${CMAKE_SOURCE_DIR})

SET_TARGET_PROPERTIES(test_ORd PROPERTIES LINKER_LANGUAGE CXX)
