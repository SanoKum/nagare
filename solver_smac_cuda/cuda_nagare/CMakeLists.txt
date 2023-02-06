find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_LIBRARY(cuda_nagare STATIC
		cudaConfig.cuh
        cudaWrapper.cu
        cudaWrapper.cuh
        calcStructualVariables_d.cu
        calcStructualVariables_d.cuh
        boundaryCond_d.cu
        boundaryCond_d.cuh
    )
	#target_link_libraries( amgcl_cuda cusparse )
	#link_directories(/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/11.8/targets/x86_64-linux/lib)
endif()
