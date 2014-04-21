include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_executable(MeshMoverDemo.out 
            MeshMoverDemo.cpp
	    )

target_link_libraries(MeshMoverDemo.out Output)                
