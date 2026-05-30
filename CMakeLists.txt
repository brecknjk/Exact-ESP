cmake_minimum_required(VERSION 3.16)
project(ShipRouting)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mbmi2")

# set(PROTOBUF_INCLUDE_DIR "/home/brecknjk/src/protobuf/src")
# set(Protobuf_PROTOC_EXECUTABLE "/home/brecknjk/src/protobuf/protoc")
set(CMAKE_BUILD_TYPE Release)
set(CGAL_USE_EPECK true)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  ) # OPTIONAL_COMPONENTS Qt5 
if ( NOT CGAL_FOUND )
  SET(CGAL_DIR "~/src/cgal")
  message(STATUS "couldn't find CGAL, trying to find it in ${CGAL_DIR}")
  find_package( CGAL QUIET COMPONENTS  ) # OPTIONAL_COMPONENTS Qt5 
  if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()  
  endif()
endif()

if ( CGAL_USE_EPECK )
  message(STATUS "Enable exact construction kernel, this will impact performance but will ensure that the results are exact.")
  add_compile_definitions(CGAL_USE_EPECK)
endif()

add_subdirectory(cpp_program src)
add_subdirectory(scripts scripts)