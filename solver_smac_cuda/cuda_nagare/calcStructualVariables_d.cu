#include "cuda_nagare/calcStructualVariables_d.cuh"

#include "flowFormat.hpp"
#include "iostream"

__global__ 
void calcStructualVariables_d 
( 
 geom_int nNormalPlanes,
 geom_int* plane_cells ,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz ,  
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz ,  
 geom_float* ccx ,  geom_float* ccy ,  geom_float* ccz ,  
 geom_float* fx   
)
{
    geom_int  ip  = blockDim.x*blockIdx.x + threadIdx.x;


    //printf("ip=%d , sx=%f , sy=%f , sz=%f \n", ip , sx[ip] , sy[ip] , sz[ip]);
    printf("sx=%d \n", sx[ip] );


    __syncthreads();

    if (ip < nNormalPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float ccx0 = ccx[ic0];
        geom_float ccy0 = ccy[ic0];
        geom_float ccz0 = ccz[ic0];

        geom_float ccx1 = ccx[ic1];
        geom_float ccy1 = ccy[ic1];
        geom_float ccz1 = ccz[ic1];

        geom_float dcx = ccx1 - ccx0;
        geom_float dcy = ccy1 - ccy0;
        geom_float dcz = ccy1 - ccy0;
        geom_float dc  = sqrtf( powf(dcx, 2.0) + pow(dcy, 2.0) + pow(dcz, 2.0));

        geom_float dc2px = pcx[ip] - ccx1;
        geom_float dc2py = pcy[ip] - ccy1;
        geom_float dc2pz = pcy[ip] - ccy1;
        geom_float dc2p  = sqrtf( powf(dc2px, 2.0) + pow(dc2py, 2.0) + pow(dc2pz, 2.0));

        fx[ip] = dc2p/dc;
    } else {
        fx[ip] = 0.0;
    }
};

void calcStructualVariables_d_wrapper(cudaConfig& cuda_cfg , mesh& msh,  variables& v)
{
    std::cout << "calcStructualVariables_d_wrapper called" << std::endl;
    calcStructualVariables_d<<<cuda_cfg.dimGrid_nplane , cuda_cfg.dimBlock>>>(
        msh.nNormalPlanes,
        msh.map_plane_cells_d,
        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"],
        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
        v.p_d["fx"] 
    );
};

