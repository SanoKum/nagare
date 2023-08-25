cmake_minimum_required(VERSION 3.18)

project(nagare CUDA CXX)

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    message("Cuda version: " ${CUDA_VERSION})
endif()

set(CMAKE_CXX_COMPILER /usr/bin/g++ )
#add_definitions( -std=c++11 -Wall )
add_definitions( -std=c++2a )
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

#cuda
set(libs_all_cuda cusparse cusolver rt dl)
set(dyn_libs_cuda amgxsh rt dl cudart_static cublas)
set(CUDA_LIB_DIR /usr/local/cuda-12.2/lib64)

#amgx
set(AMGX_INCLUDES ${THRUST_DIR} /home/kumpei/src/AMGX/include )
set(AMGX_LIB_DIR /home/kumpei/src/AMGX/build)

#hdf
set(hdf5_inc /home/kumpei/src/HighFive/include /home/kumpei/src/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/include)
set(hdf5_libdir /home/kumpei/src/hdf/HDF5-1.14.1-Linux/HDF_Group/HDF5/1.14.1/lib)

include_directories(./ ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${hdf5_inc} ${AMGX_INCLUDES})
link_directories(./ ${hdf5_libdir} ${AMGX_LIB_DIR} ${CUDA_LIB_DIR})
#link_directories(./ ${hdf5_libdir} ${AMGX_LIB_DIR} )

# --------------
# *** nagare ***
# --------------
add_executable(nagare main.cpp gradient.cpp mesh/mesh.cpp 
               output/output.cpp update.cpp variables.cpp 
               boundaryCond.cpp solvePoisson_amgx.cpp solveNavierStokes.cpp 
               input/solverConfig.cpp calcCFL.cpp setStructualVariables.cpp)
#target_include_directories(nagare PUBLIC )
target_link_libraries(nagare hdf5 yaml-cpp ${dyn_libs_cuda})
#target_link_libraries(nagare hdf5 yaml-cpp ${libs_all_cuda})

# ----------------------
# *** mesh conversion***
# ----------------------
add_executable(convertGmshToNagare mesh/convertGmshToNagare.cpp mesh/mesh.cpp variables.cpp boundaryCond.cpp input/solverConfig.cpp )
target_link_libraries(convertGmshToNagare hdf5 yaml-cpp)
