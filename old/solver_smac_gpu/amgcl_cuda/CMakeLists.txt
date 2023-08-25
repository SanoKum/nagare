find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
    CUDA_ADD_LIBRARY(amgcl_cuda STATIC
        solvePoisson_amgcl_cuda.cu
        solvePoisson_amgcl_cuda.h
    )
	target_link_libraries( amgcl_cuda ${CUDA_cusparse_LIBRARY} )
	#target_link_libraries( amgcl_cuda cusparse )
	#link_directories(/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/11.8/targets/x86_64-linux/lib)
endif()
