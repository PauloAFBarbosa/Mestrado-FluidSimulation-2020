file(GLOB sources src/*.cpp src/*.cu)
file(GLOB headers src/*.h)

add_library(nauPassCudaHash SHARED
	${sources} ${headers})
	
add_dependencies(nauPassCudaHash nau)
	
include_directories( 
	../../nau/src
	../../nau/src/nau
)	

if (WIN32)
	add_definitions(-D_WINDLL)
endif(WIN32)
	
target_link_libraries( nauPassCudaHash $<TARGET_LINKER_FILE:nau> glbinding ${CUDA_LIBRARIES})
include_directories(${CUDA_INCLUDE_DIRS})


install (TARGETS nauPassCudaHash DESTINATION lib)
install(FILES ${headers} DESTINATION include)	

add_custom_command(
	TARGET nauPassCudaHash POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
	$<TARGET_FILE:nauPassCudaHash>
    ${CMAKE_BINARY_DIR}/bin/nauSettings/plugins/$<CONFIGURATION>/pass/$<TARGET_FILE_NAME:nauPassCudaHash>)
					   
if (WIN32)
	set_target_properties (nauPassCudaHash 
		PROPERTIES FOLDER plugins)	
endif(WIN32)		
