project(shrinker)

set(LINK_LIBS
  ${LINK_LIBS}
)

set(HEADERS
   ${HEADERS}
)
set(SOURCES_CPP
	${SOURCES_CPP}
	${CMAKE_CURRENT_SOURCE_DIR}/shrinkerNew.cpp
)

SET(CGAL_DIR "~/src/cgal")
# CGAL and its components
find_package( CGAL QUIET COMPONENTS  ) # OPTIONAL_COMPONENTS Qt5

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
# else()
#   add_definitions(-DUSE)
endif()

add_executable(${PROJECT_NAME} ${SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if ( CGAL_FOUND )
  target_link_libraries(shrinker CGAL::CGAL)
else()
    message(STATUS "couldn't link CGAL")
endif()