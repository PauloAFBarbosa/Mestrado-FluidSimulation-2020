file(GLOB sources src/*.cpp)
file(GLOB headers src/*.h)

add_library(nauPassCuda SHARED
	${sources} ${headers})
	
add_dependencies(nauPassCuda nau)
	
include_directories( 
	../../nau/src
	../../nau/src/nau
)	

if (WIN32)
	add_definitions(-D_WINDLL)
endif(WIN32)
	
target_link_libraries( nauPassCuda $<TARGET_LINKER_FILE:nau> glbinding ${CUDA_LIBRARIES})
include_directories(${CUDA_INCLUDE_DIRS})


install (TARGETS nauPassCuda DESTINATION lib)
install(FILES ${headers} DESTINATION include)	

add_custom_command(
	TARGET nauPassCuda POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
	$<TARGET_FILE:nauPassCuda>
    ${CMAKE_BINARY_DIR}/bin/nauSettings/plugins/$<CONFIGURATION>/pass/$<TARGET_FILE_NAME:nauPassCuda>)
					   
if (WIN32)
	set_target_properties (nauPassCuda 
		PROPERTIES FOLDER plugins)	
endif(WIN32)		
