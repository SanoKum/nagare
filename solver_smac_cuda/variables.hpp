#pragma once
#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "cuda_nagare/cudaConfig.cuh"

class variables {
public:
    std::map<std::string, std::vector<flow_float>> c;
    std::map<std::string, std::vector<flow_float>> p;
    std::map<std::string, flow_float*> c_d;
    std::map<std::string, flow_float*> p_d;

    const std::list<std::string> cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , "deltaP" ,
        "roN" , "UxN" , "UyN" , "UzN" , "TN" ,
        "dTdx", "dTdy", "dTdz",
        "dUxdx", "dUxdy", "dUxdz",
        "dUydx", "dUydy", "dUydz",
        "dUzdx", "dUzdy", "dUzdz",
        "dPdx", "dPdy", "dPdz",
        "divU*vol" , "divU" , "divU_star",
        "convx" , "convy" , "convz" , "convT",
        "diffx" , "diffy" , "diffz" , "diffT",
        "dP"    , "dPPdx", "dPPdy", "dPPdz",
        "cfl"   ,

        // Mesh Structure
        "volume" , "ccx" , "ccy" , "ccz"

    };

    const std::list<std::string> planeValNames = 
    {
        "US" , "T" ,  
        "USN", "Ux" , "Uy" , "Uz", "ro" , "P",

        // Mesh Structure
        "sx" , "sy" , "sz" , "ss" ,
        "pcx" , "pcy" , "pcz" ,
        "fx" , 
    };

    const std::list<std::string> output_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , 
        "divU" , "divU_star",
        "convx" , "convy" , "convz" , 
        "diffx" , "diffy" , "diffz" , 
        "dPPdx" , "dPPdy"   , "dPPdz",
        "dP"  ,
        "cfl"
    };

    //not yet implemented
    const std::list<std::string> read_cellValNames = 
    {
        "ro"  , "Ux"  , "Uy"  , "Uz"  , "T" , "P" , "dP"
    };

    variables();

    variables(mesh&);
    ~variables();

    void copyVariables_cell_plane_H2D();

    void copyVariables_cell_plane_D2H();

    void setStructualVariables_d(cudaConfig& cuda_cfg , mesh& msh , variables& v);
};
