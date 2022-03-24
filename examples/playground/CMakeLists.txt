SET(TESTNAMEV playground1)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "playground1")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV playground2)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main2.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "playground2")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV playground3)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main3.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "playground3")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})


SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV playground4)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} main4.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "playground4")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)


SET(TESTNAMEV gmsh2exo)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} gmsh2exo.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "gmsh2exo")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)


SET(TESTNAMEV gmsh2exo_LA)
message("=== Adding Example : ${TESTNAMEV}")
add_executable(${TESTNAMEV} gmsh2exo_LA.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "gmsh2exo_LA")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV transfer_fibers)
message("=== Adding Example : ${TESTNAMEV}")

add_executable(${TESTNAMEV} transfer_fibers.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "transfer_fibers")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV remove_valves)
message("=== Adding Example : ${TESTNAMEV}")
add_executable(${TESTNAMEV} remove_valves.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "remove_valves")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)



SET(TESTNAMEV benchmark_fibers)
message("=== Adding Example : ${TESTNAMEV}")
add_executable(${TESTNAMEV} benchmark_fibers.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "remove_valves")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)

SET(TESTNAMEV read_exodus_fibers)
message("=== Adding Example : ${TESTNAMEV}")
add_executable(${TESTNAMEV} read_exodus_fibers.cpp)

set_target_properties(${TESTNAMEV} PROPERTIES  OUTPUT "read_exodus_fibers")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${LIBMESH_DIR}/include")

target_link_libraries(${TESTNAMEV} beatit)
target_link_libraries(${TESTNAMEV} ${LIBMESH_LIB})
target_link_libraries(${TESTNAMEV} ${VTK_LIBRARY_DIRS})

SET_TARGET_PROPERTIES(${TESTNAME} PROPERTIES LINKER_LANGUAGE CXX)
