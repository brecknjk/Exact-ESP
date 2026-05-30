project(extractor)

if ( USE_OSMPBF )
  set(LINK_LIBS
    ${LINK_LIBS}
    osmpbf
  )
  add_definitions(-DUSE_OSMPBF)
else()
  set(LINK_LIBS
	  ${LINK_LIBS}
  )
endif()

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/osmExtractor.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/triangulator.hpp
)
set(SOURCES_CPP
	${SOURCES_CPP}
	${CMAKE_CURRENT_SOURCE_DIR}/extractor.cpp
)

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

add_executable(${PROJECT_NAME} ${SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if ( CGAL_FOUND )
  target_link_libraries(extractor CGAL::CGAL)
else()
    message(STATUS "couldn't link CGAL")
endif()