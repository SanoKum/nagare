find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_LIBRARY(thrust_nagare STATIC
        variables_thrust.cu
        variables_thrust.h
        setInitial_thrust.cu
        setInitial_thrust.h
        update_thrust.cu
        update_thrust.h
        copyDeviceToHost_cell.cu
        copyDeviceToHost_cell.h
    )
	#target_link_libraries( amgcl_cuda cusparse )
	#link_directories(/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/11.8/targets/x86_64-linux/lib)
endif()
