file(GLOB sources src/*.cpp)
file(GLOB headers src/*.h)

add_library(nauSortPlugin SHARED
	${sources} ${headers})
	
add_dependencies(nauSortPlugin nau)
	
include_directories( 
	../../nau/src
	../../nau/src/nau
)	

if (WIN32)
	add_definitions(-D_WINDLL)
endif(WIN32)
	
target_link_libraries( nauSortPlugin $<TARGET_LINKER_FILE:nau> glbinding)
	
install (TARGETS nauSortPlugin DESTINATION lib)
install(FILES ${headers} DESTINATION include)	

add_custom_command(
	TARGET nauSortPlugin POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
	$<TARGET_FILE:nauSortPlugin>
    ${CMAKE_BINARY_DIR}/bin/nauSettings/plugins/$<CONFIGURATION>/pass/$<TARGET_FILE_NAME:nauSortPlugin>)
					   
if (WIN32)
	set_target_properties (nauSortPlugin 
		PROPERTIES FOLDER plugins)	
endif(WIN32)		
