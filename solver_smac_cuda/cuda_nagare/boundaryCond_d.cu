#include "cuda_nagare/boundaryCond_d.cuh"

__global__ 
void slip_d 
( 
 geom_int ibSt, geom_int ibEd,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P
)
{
    geom_int ibl  = blockDim.x*blockIdx.x + threadIdx.x;
    geom_int ib = ibSt + ibl;

    if (ib < ibEd) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        geom_float Un =  (sx[ip]*Uxp[ic] + sy[ip]*Uyp[ic] + sz[ip]*Uzp[ic])/ss[ip];

        Tp[ip]  = T[ic];
        Uxp[ip] = Ux[ic] - Un*sx[ip]/ss[ip];
        Uyp[ip] = Uy[ic] - Un*sy[ip]/ss[ip];
        Uzp[ip] = Uz[ic] - Un*sz[ip]/ss[ip];
        rop[ip] = ro[ic];
        Pp[ip]  = P[ic];
    }
};

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    std::cout << "aaa" << std::endl;
    std::cout << "T_d[0] = " << var.p_d["T"][0] << std::endl;

    slip_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        bc.iBPlanes[0] , bc.iBPlanes.back(), 
        msh.map_bplane_plane_d,  msh.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["T"]  , var.c_d["T"] ,
        var.p_d["Ux"] , var.c_d["Ux"] ,
        var.p_d["Uy"] , var.c_d["Uy"] ,
        var.p_d["Uz"] , var.c_d["Uz"] ,
        var.p_d["ro"] , var.c_d["ro"] ,
        var.p_d["P"]  , var.c_d["P"] 
    );
}

__global__ 
void wall_d 
( 
 geom_int ibSt, geom_int ibEd,
 // mesh structure
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* Tp  ,  flow_float* T   ,
 flow_float* Uxp ,  flow_float* Ux  ,
 flow_float* Uyp ,  flow_float* Uy  ,
 flow_float* Uzp ,  flow_float* Uz  ,
 flow_float* rop ,  flow_float* ro  ,
 flow_float* Pp  ,  flow_float* P   ,
 // bvar
 flow_float* Uxb ,  
 flow_float* Uyb ,  
 flow_float* Uzb   
)
{
    geom_int ibl  = blockDim.x*blockIdx.x + threadIdx.x;
    geom_int ib = ibSt + ibl;

    if (ib < ibEd) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];

        Uxp[ip] = Uxb[ib];
        Uyp[ip] = Uyb[ib];
        Uzp[ip] = Uzb[ib];
        
        Tp[ip]  = T[ic];
        rop[ip] = ro[ic];
        Pp[ip]  = P[ic];
    }
};

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        bc.iBPlanes[0] , bc.iBPlanes.back(), 
        msh.map_bplane_plane_d,  msh.map_bplane_cell_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.p_d["Tp"]  , var.c_d["T"] ,
        var.p_d["Uxp"] , var.c_d["Ux"] ,
        var.p_d["Uyp"] , var.c_d["Uy"] ,
        var.p_d["Uzp"] , var.c_d["Uz"] ,
        var.p_d["rop"] , var.c_d["ro"] ,
        var.p_d["Pp"]  , var.c_d["P"]  , 

        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"]
    );
}



