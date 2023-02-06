#include <iostream>
#include <vector>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "output/output.hpp"
#include "boundaryCond.hpp"

#include "setStructualVariables.hpp"

#include "gradient.hpp"

#include "solvePoisson_amgcl.hpp"
#include "solveNavierStokes.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

#include "calcCFL.hpp"

// cuda
#include "cuda_nagare/cudaWrapper.cuh"
#include "cuda_nagare/cudaConfig.cuh"

int main(void) {
    clock_t start = clock();

    cout << "Read Solver Config \n";
    solverConfig cfg; 
    cfg.read("solverConfig.yaml");

    cout << "Read Mesh \n";
    mesh msh;
    if (cfg.meshFormat == "hdf5") {
        msh.readMesh(cfg.meshFileName);
    } else if (cfg.meshFormat == "gmsh") {
        gmshReader gmsh = gmshReader(cfg.meshFileName);
        msh = gmsh.getMesh();
    } else {
        cerr << "Error unknown mesh format: " << cfg.meshFormat << endl;
        return 1;
    }

    cout << "Init Matrix \n";
    matrix mat_poi = matrix();
    mat_poi.initMatrix(msh);

    matrix mat_ns = matrix();
    mat_ns.initMatrix(msh);

    cout << "Read Boundary Conditions \n";
    readAndSetBcondConfig(msh.bconds);

    cout << "Set Initial Values \n";
    variables v = variables(msh);
    setInitial(cfg , msh , v);

    updateVariablesForNextLoop(cfg , msh , v , mat_poi);
    v.copyVariables_cell_plane_H2D();

    // prepare for cuda
    cudaConfig cuda_cfg = cudaConfig(msh);
    msh.setMeshMap_d();

    v.setStructualVariables_d( cuda_cfg , msh , v);

    cout << "before Tp0=" << v.p["T"][0] << endl;
    cout << "before Tc0=" << v.c["T"][0] << endl;
    cout << "before Tc2=" << v.c["T"][2] << endl;
    cout << "before sx0=" << v.p["sx"][0] << endl;
    cout << "before sy0=" << v.p["sy"][0] << endl;
    cout << "before sz0=" << v.p["sz"][0] << endl;


    v.copyVariables_cell_plane_D2H();

    cout << "Tp0=" << v.p["T"][0] << endl;
    cout << "Tc0=" << v.c["T"][0] << endl;
    cout << "Tc2=" << v.c["T"][2] << endl;
    cout << "sx0=" << v.p["sx"][0] << endl;
    cout << "sy0=" << v.p["sy"][0] << endl;
    cout << "sz0=" << v.p["sz"][0] << endl;

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

//        cout << " Step =" << iStep << "\n";
//        setBcondsValue(cfg , cuda_cfg , msh , v , mat_poi);
//temp
//temp        calcGradient(cfg , msh , v);
//temp
//temp        solveNavierStokes(cfg , msh , v , mat_ns);
//temp
//temp        setMatrixPoisson(cfg , msh , v , mat_poi);
//temp
//temp        if (cfg.gpu == 0) {
//temp            solvePoisson(cfg , msh , v , mat_poi);
//temp        } else if (cfg.gpu == 1) {
//temp            callAmgclCuda(cfg, msh, v, mat_poi );
//temp        } else {
//temp            cerr << "unknown gpu option : " << cfg.gpu << endl;
//temp            return 1;//tem
//temp        }
//temp
//temp        correctPresVel(cfg , msh , v);
//temp
//temp        // *** output hdf5 ***
//temp        outputH5_XDMF(cfg , msh, v, iStep);
//temp
//temp        // *** update TN, diffT, etc.***
//temp        updateVariablesForNextLoop(cfg , msh , v , mat_poi);
//temp
//temp        calcCFL(cfg , msh , v);
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 

	return 0;
}