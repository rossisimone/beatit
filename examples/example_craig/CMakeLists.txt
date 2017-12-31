SET(TESTNAMEV example_craig)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_craig")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data.beat  IS_NEWER_THAN ${GetPotFile} )


SET(TESTNAMEV compute_activation_times)
message("=== Adding Example : ${TESTNAMEV}")


add_executable(${TESTNAMEV} main_compute_activation_times.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "compute_activation_times")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data_at.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data_at.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data_at.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data_at.beat  IS_NEWER_THAN ${GetPotFile} )


SET(TESTNAMEV compute_fibers_pw)
message("=== Adding Example : ${TESTNAMEV}")


add_executable(${TESTNAMEV} main_fibers.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "compute_fibers_pw")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET(GetPotFile "${CMAKE_CURRENT_BINARY_DIR}/data_fibers.beat")
IF ( ${CMAKE_CURRENT_SOURCE_DIR}/data_fibers.beat  IS_NEWER_THAN ${GetPotFile} )
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/data_fibers.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
ENDIF (${CMAKE_CURRENT_SOURCE_DIR}/data_fibers.beat  IS_NEWER_THAN ${GetPotFile} )