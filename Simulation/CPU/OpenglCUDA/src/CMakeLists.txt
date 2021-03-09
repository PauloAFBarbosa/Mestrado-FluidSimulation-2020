cmake_minimum_required(VERSION 3.6)

# Project Name
PROJECT(Octree)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_executable(${PROJECT_NAME} main.cpp )

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})

if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)



if  (WIN32)

	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Octree)
	
	message(STATUS "GLEW_DIR set to: " ${GLEW_DIR})
	set(GLEW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../toolkits/glew" CACHE PATH "Path to glew folder")
	if (NOT EXISTS "${GLEW_DIR}/GL/glew.h" OR NOT EXISTS "${GLEW_DIR}/glew32.lib")
		message(ERROR ": glew not found")
	endif (NOT EXISTS "${GLEW_DIR}/GL/glew.h" OR NOT EXISTS "${GLEW_DIR}/glew32.lib")	
	include_directories(${GLEW_DIR})
	link_directories(${GLEW_DIR})

	message(STATUS "GLUT_DIR set to: " ${GLUT_DIR})
	set(GLUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../toolkits/glut" CACHE PATH "Path to GLUT folder")
	if (NOT EXISTS "${GLUT_DIR}/GL/glut.h" OR NOT EXISTS "${GLUT_DIR}/glut32.lib")
		message(ERROR ": GLUT not found")
	endif (NOT EXISTS "${GLUT_DIR}/GL/glut.h" OR NOT EXISTS "${GLUT_DIR}/glut32.lib")	
	include_directories(${GLUT_DIR})
	link_directories(${GLUT_DIR})

	

	

	target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLEW_DIR}/glew32.lib ${GLUT_DIR}/glut32.lib )
	
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