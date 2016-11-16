
#set(HEADERS 
#    ${HEADERS} 
#    "${CMAKE_CURRENT_SOURCE_DIR}/io/*.hpp")

#set(SOURCES 
#    ${SOURCES} 
#    "${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp")


file(GLOB_RECURSE SOURCES *.cpp)
file(GLOB_RECURSE HEADERS  *.hpp)

#add_subdirectory(io)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library( beatit  ${HEADERS} ${SOURCES})
target_link_libraries (beatit ${LIBMESH_LIB})