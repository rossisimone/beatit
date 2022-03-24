SET(TESTNAME example_josue)
message("=== Adding Example : ${TESTNAME}")

add_executable(${TESTNAME} main.cpp)

set_target_properties(${TESTNAME} PROPERTIES  OUTPUT "example_josue")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")


target_link_libraries(${TESTNAME} beatit)
target_link_libraries(${TESTNAME} ${LIBMESH_LIB})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/em.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/electrophysiology.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/activation.beat  ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

SET(TESTNAMEV example_josue_bath)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main_bath.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_josue_bath")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV example_josue_bath_convergence)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main_bath_convergence.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "example_josue_bath_convergence")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)