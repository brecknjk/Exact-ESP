project(router)

set(LINK_LIBS
	${LINK_LIBS}
)
set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/edgeLocator.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/polyanya.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/deadendFinder.hpp
   ${CMAKE_CURRENT_SOURCE_DIR}/graph.hpp
)
set(SOURCES_CPP
	${SOURCES_CPP}
	${CMAKE_CURRENT_SOURCE_DIR}/router.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/graph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deadendFinder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/edgeLocator.cpp
)

SET(CGAL_DIR "~/src/cgal")
# CGAL and its components
find_package( CGAL QUIET COMPONENTS  ) # OPTIONAL_COMPONENTS Qt5 

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

create_single_source_cgal_program( "router.cpp" )

# add_executable(${PROJECT_NAME} ${SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if ( OpenMP_CXX_FOUND )
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# if ( CGAL_Qt5_FOUND )
#   target_link_libraries(${PROJECT_NAME} PUBLIC CGAL::CGAL_Basic_viewer)
# endif()