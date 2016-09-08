add_executable(test_ORd main.cpp)

set_target_properties(test_ORd PROPERTIES  OUTPUT "test_ORd")


target_link_libraries(test_ORd beatit)
target_link_libraries(test_ORd ${LIBMESH_LIB})

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

SET_TARGET_PROPERTIES(test_ORd PROPERTIES LINKER_LANGUAGE CXX)

SET(OctaveFile "${CMAKE_CURRENT_BINARY_DIR}/plot_variables.m")

IF ( ${CMAKE_CURRENT_SOURCE_DIR}/plot_variables.m  IS_NEWER_THAN ${OctaveFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/plot_variables.m  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/plot_variables.m  IS_NEWER_THAN ${OctaveFile} )
    
     
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_variables.m DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


SET(TESTNAME test_beatit_ORd)
add_executable(${TESTNAME} main_ORd.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "test_beatit_ORd")

target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})
include_directories ("${PROJECT_SOURCE_DIR}/src")

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

add_test(${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/test_beatit_ORd)
