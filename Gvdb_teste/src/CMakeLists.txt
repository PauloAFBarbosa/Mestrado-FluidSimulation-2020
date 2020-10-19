cmake_minimum_required(VERSION 3.6)

# Project Name
PROJECT(Gvdb_teste)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_executable(${PROJECT_NAME} main.cpp Octree.h Octree.cpp Point.h Point.cpp)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})

if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)



if  (WIN32)

	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Gvdb_teste)
	message(STATUS "GLUT_DIR set to: " ${GLUT_DIR})
	set(GLUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../toolkits/glut" CACHE PATH "Path to GLUT folder")
	if (NOT EXISTS "${GLUT_DIR}/GL/glut.h" OR NOT EXISTS "${GLUT_DIR}/glut32.lib")
		message(ERROR ": GLUT not found")
	endif (NOT EXISTS "${GLUT_DIR}/GL/glut.h" OR NOT EXISTS "${GLUT_DIR}/glut32.lib")	
	include_directories(${GLUT_DIR})
	link_directories(${GLUT_DIR})

	# If the gvdb target hasn't yet been defined (for instance, if this is
	# being called from a build_all script), import the GVDB project directly:
	if(NOT TARGET gvdb)
	    # One fast way to add GVDB as a build dependency of the project is to add its subdirectory as follows:
	    add_subdirectory("gvdb_library/" "$CMAKE_CURRENT_BINARY_DIR}/gvdb_library")
	    # Alternatively, if GVDB has already been built, then we can import its GVDBConfig.cmake list file:
	    # find_package(GVDB CONFIG REQUIRED
	    #    PATHS "${CMAKE_CURRENT_BINARY_DIR}/../gvdb_library")
	endif()

	file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../toolkits/gvdb/" DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/../Build/")

	target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_DIR}/glut32.lib gvdb )
	
else (WIN32) #Linux and Mac

	find_package(GLUT REQUIRED)
	include_directories(${GLUT_INCLUDE_DIR})
	link_directories(${GLUT_LIBRARY_DIRS})
	add_definitions(${GLUT_DEFINITIONS})
	target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
	if(NOT GLUT_FOUND)
	   message(ERROR ": GLUT not found!")
	endif(NOT GLUT_FOUND)
	
endif(WIN32)