# minimum required version of CMake
cmake_minimum_required(VERSION 2.8.9)

# Project name
project(BeatIt)
set (BeatIt_VERSION_MAJOR 0)
set (BeatIt_VERSION_MINOR 0)




# Include source directories
add_subdirectory(src)
add_subdirectory(testsuite)

#include_directories(src)
include_directories ("${PROJECT_SOURCE_DIR}/src")

# set flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# 

# Packages 
find_package(MPI REQUIRED)

#
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
message(" Module path: ${CMAKE_MODULE_PATH}")
find_package(LIBMESH REQUIRED)


# Summary 
message(" Compiling in ${CMAKE_BUILD_TYPE}")


get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()
