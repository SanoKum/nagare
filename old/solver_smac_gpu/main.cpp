#include <iostream>
#include <vector>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <sstream>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"
#include "output/output.hpp"

#include "mesh/gmshReader.hpp"
#include "setStructualVariables.hpp"
#include "boundaryCond.hpp"

#include "gradient.hpp"

#include "solvePoisson_amgcl.hpp"
#include "solveNavierStokes.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

#include "calcCFL.hpp"

// thrust
#include "thrust_nagare/variables_thrust.h"
#include "thrust_nagare/setInitial_thrust.h"
#include "thrust_nagare/update_thrust.h"
#include "thrust_nagare/copyDeviceToHost_cell.h"
#include "thrust_nagare/setStructualVariables_thrust.h"

int main(void) {
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
    variables_d v_d = variables_d(msh, v.cellValNames , v.planeValNames); // device

    setInitial(cfg , msh , v);
    setInitial_thrust(cfg , msh , v_d);

    setStructualVariables(cfg , msh , v);
    setStructualVariables_thrust(cfg , msh , v_d);

    updateVariablesForNextLoop(cfg , msh , v , mat_poi);
    updateVariablesForNextLoop_thrust(cfg , msh , v_d , mat_poi);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep <cfg.nStep ; iStep++) {

        cout << " Step =" << iStep << "\n";
        setBcondsValue(cfg , msh , v , mat_poi);

        calcGradient(cfg , msh , v);

        solveNavierStokes(cfg , msh , v , mat_ns);

        setMatrixPoisson(cfg , msh , v , mat_poi);

        if (cfg.gpu == 0) {
            solvePoisson(cfg , msh , v , mat_poi);
        } else if (cfg.gpu == 1) {
            callAmgclCuda(cfg, msh, v, mat_poi );
        } else {
            cerr << "unknown gpu option : " << cfg.gpu << endl;
            return 1;
        }

        correctPresVel(cfg , msh , v);


        // *** output hdf5 ***
        copyDeviceToHost_cell(v_d , v);
        outputH5_XDMF(cfg , msh, v, iStep);

        // *** update TN, diffT, etc.***
        updateVariablesForNextLoop(cfg , msh , v , mat_poi);
        updateVariablesForNextLoop_thrust(cfg , msh , v_d , mat_poi);

        calcCFL(cfg , msh , v);
    }
    
	return 0;
}