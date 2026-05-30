project(utils)

set(SOURCE
   ${SOURCE}
)
set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/mathLib.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/io.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/types.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/settings.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/section.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/timer.hpp
)
set(SOURCES_CPP
	${SOURCES_CPP}
	${CMAKE_CURRENT_SOURCE_DIR}/route_to_gl.cpp
)

SET(CGAL_DIR "~/src/cgal")
# CGAL and its components
find_package( CGAL QUIET COMPONENTS  ) # OPTIONAL_COMPONENTS Qt5 

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

add_executable(route_to_gl ${SOURCES_CPP})
target_link_libraries(route_to_gl ${LINK_LIBS})
target_include_directories(route_to_gl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if ( CGAL_FOUND )
  target_link_libraries(route_to_gl CGAL::CGAL)
else()
    message(STATUS "couldn't link CGAL")
endif()