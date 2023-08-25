#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_nagare/cudaWrapper.cuh"
#include "cuda_nagare/calcStructualVariables_d.cuh"

variables::variables() {};

variables::variables(mesh& msh)
{
    for (auto& cellValName : cellValNames)
    {
        this->c[cellValName].resize(msh.nCells);

        //cudaWrapper::cudaMalloc_wrapper(&(this->c_d[cellValName]) , msh.nCells);
        CHECK_CUDA_ERROR( cudaMalloc((void**) &(this->c_d[cellValName]), msh.nCells*sizeof(flow_float)) );

    }
    for (auto& planeValName : planeValNames)
    {
        this->p[planeValName].resize(msh.nPlanes);

        //cudaWrapper::cudaMalloc_wrapper(&(this->p_d[planeValName]) , msh.nPlanes);
        CHECK_CUDA_ERROR( cudaMalloc((void**) &(this->p_d[planeValName]), msh.nPlanes*sizeof(flow_float)) );
    }
}

variables::~variables()
{

    for (auto& cellValName : cellValNames)
    {
        CHECK_CUDA_ERROR( cudaFree(this->c_d[cellValName]) );
    }

    for (auto& planeValName : planeValNames)
    {
        CHECK_CUDA_ERROR( cudaFree(this->p_d[planeValName]) );
    }

}

void variables::copyVariables_cell_plane_H2D_all()
{
    for (auto& name : this->cellValNames)
    {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->c[name].data() , this->c_d[name], this->c[name].size());
    }
    for (auto& name : this->planeValNames)
    {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->p[name].data() , this->p_d[name], this->p[name].size());
    }
}

void variables::copyVariables_cell_H2D(std::string name)
{
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->c[name].data() , this->c_d[name], this->c[name].size());
}
void variables::copyVariables_plane_H2D(std::string name)
{
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->p[name].data() , this->p_d[name], this->c[name].size());
}

void variables::copyVariables_cell_plane_D2H_all()
{
    for (auto& name : this->cellValNames)
    {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d[name], this->c[name].data() , this->c[name].size());
    }
    for (auto& name : this->planeValNames)
    {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d[name], this->p[name].data() , this->p[name].size());
    }
}

void variables::copyVariables_cell_D2H(std::string name)
{
    cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d[name], this->c[name].data(), this->c[name].size());
}

void variables::copyVariables_plane_D2H(std::string name)
{
    cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d[name], this->p[name].data(), this->p[name].size());
}

void variables::setStructualVariables_d(cudaConfig& cuda_cfg , mesh& msh )
{
    geom_float* sx;
    geom_float* sy;
    geom_float* sz;
    geom_float* ss;
    geom_float* pcx;
    geom_float* pcy;
    geom_float* pcz;
    geom_float* ccx;
    geom_float* ccy;
    geom_float* ccz;
    geom_float* volume;

    sx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    ss = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);

    ccx = (geom_float*)malloc(sizeof(geom_float)*msh.nCells);
    ccy = (geom_float*)malloc(sizeof(geom_float)*msh.nCells);
    ccz = (geom_float*)malloc(sizeof(geom_float)*msh.nCells);
    volume = (geom_float*)malloc(sizeof(geom_float)*msh.nCells);


    for (geom_int ip=0; ip<msh.nPlanes; ip++)
    {
        sx[ip] = msh.planes[ip].surfVect[0];
        sy[ip] = msh.planes[ip].surfVect[1];
        sz[ip] = msh.planes[ip].surfVect[2];
        ss[ip] = msh.planes[ip].surfArea;
        pcx[ip] = msh.planes[ip].centCoords[0];
        pcy[ip] = msh.planes[ip].centCoords[1];
        pcz[ip] = msh.planes[ip].centCoords[2];
    }

    for (geom_int ic=0; ic<msh.nCells; ic++)
    {
        ccx[ic] = msh.cells[ic].centCoords[0];
        ccy[ic] = msh.cells[ic].centCoords[1];
        ccz[ic] = msh.cells[ic].centCoords[2];
        volume[ic] = msh.cells[ic].volume;
    }

    gpuErrchk(cudaMemcpy(this->p_d["sx"] , sx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice));
    cudaMemcpy(this->p_d["sy"] , sy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["sz"] , sz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["ss"] , ss , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(this->p_d["pcx"] , pcx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["pcy"] , pcy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->p_d["pcz"] , pcz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(this->c_d["ccx"] , ccx , msh.nCells*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["ccy"] , ccy , msh.nCells*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["ccz"] , ccz , msh.nCells*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(this->c_d["volume"] , volume , msh.nCells*sizeof(geom_float) , cudaMemcpyHostToDevice);

    calcStructualVariables_d_wrapper(cuda_cfg , msh , *this);

    //for (geom_int ip=0; ip<msh.nPlanes; ip++)
    //{
    //    printf("ip=%d , sx=%f , sy=%f , sz=%f\n", ip, sx[ip], sy[ip], sz[ip]);
    //}


    free(sx) ; free(sy) ; free(sz) ; free(ss);
    free(pcx); free(pcy); free(pcz);
    free(ccx); free(ccy); free(ccz);
    free(volume); 

}