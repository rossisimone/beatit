# minimum required version of CMake
cmake_minimum_required(VERSION 2.8.9)

# Project name
project(BeatIt)
set (BeatIt_VERSION_MAJOR 0)
set (BeatIt_VERSION_MINOR 0)

# REQUIRES C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("\n\t Compiling in ${CMAKE_BUILD_TYPE}\n")
if( ${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  message("\t Adding -D_GLIBCXX_DEBUG to the CXX flags:\n")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -Wall -D_GLIBCXX_DEBUG")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -ffast-math")  
endif()
message(STATUS "CXX flags: ${CMAKE_CXX_FLAGS} \n")


# Packages 
find_package(MPI REQUIRED)
find_package(VTK REQUIRED NO_MODULE)
message("-- VTK_DIR: ${VTK_DIR}")
message("-- VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")
message("-- VTK_LIBRARY_DIRS: ${VTK_LIBRARY_DIRS}")
message("-- VTK_LIBRARIES: ${VTK_LIBRARIES}")

find_package(Boost REQUIRED)
message("-- BOOST_LIBRARYDIR: ${BOOST_LIBRARYDIR}")
message("-- Boost_LIBRARIES: ${Boost_LIBRARIES}")
message("-- Boost_COMPILER: ${Boost_COMPILER}")

#
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
message("-- Module path: ${CMAKE_MODULE_PATH}")
find_package(LIBMESH REQUIRED)

set(LIBMESH_LIB ${LIBMESH_LIBRARIES})




#include_directories(src)
include_directories ("${Boost_INCLUDE_DIR}")

if("${VTK_INCLUDE_DIRS}" STREQUAL "")
message("-- Not including VTK")
else("${VTK_INCLUDE_DIRS}" STREQUAL "")
message("-- Including VTK")
include_directories ("${VTK_INCLUDE_DIRS}")
endif("${VTK_INCLUDE_DIRS}" STREQUAL  "")

include_directories ("${LIBMESH_DIR}/include")
message("-- PETSC_INCLUDE_DIRS: ${PETSC_DIR}/include")
message("-- PETSC_LIBRARIES: ${PETSC_LIBRARIES}")
include_directories ("${PETSC_DIR}/include")
include_directories ("${PROJECT_SOURCE_DIR}/src")
message("-- MPI LIBRARIES: ${MPI_LIBRARIES}")
message("-- Boost Include: ${Boost_INCLUDE_DIR}")
message("-- Boost Libraries: ${Boost_LIBRARY_DIR}")

message("-- EIGEN Directory: ${EIGEN_DIR}")
include_directories ("${EIGEN_DIR}")

enable_testing()

# Include source directories
message("-- Including src folder")
add_subdirectory(src)
message("-- Adding library")

file(GLOB_RECURSE SOURCES src/*.cpp)
file(GLOB_RECURSE HEADERS src/*.hpp)
add_library( beatit  ${HEADERS} ${SOURCES})
target_link_libraries (beatit ${LIBMESH_LIB})
target_link_libraries (beatit ${PETSC_LIBRARIES})
target_link_libraries (beatit ${VTK_LIBRARIES})

add_subdirectory(examples)
add_subdirectory(testsuite)

# set flags
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")




# Summary 

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

message("\n\tPlease type `make` in order to compile BeatIt.\n")


message(STATUS "CXX DEBUG flags: ${CMAKE_CXX_FLAGS} \n")
