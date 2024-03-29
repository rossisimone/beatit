SET(TESTNAMEV example_san_atria)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_san_atria")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})
target_link_libraries(${TESTNAMEV} ${PETSC_LIBRARIES})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/monowave.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  IS_NEWER_THAN ${GetPotFile} )


SET(TESTNAMEV example_fibrosis)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main_fibrosis.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_fibrosis")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})
target_link_libraries(${TESTNAMEV} ${PETSC_LIBRARIES})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/monowave.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/monowave.beat  IS_NEWER_THAN ${GetPotFile} )
