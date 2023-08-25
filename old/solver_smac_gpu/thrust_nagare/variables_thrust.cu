#include "mesh/mesh.hpp"

#include "thrust_nagare/variables_thrust.h"

//#include <list>
//#include <vector>
//#include <string>

variables_d::variables_d(mesh& msh, std::list<std::string> cellValNames , std::list<std::string> planeValNames)
{
    for (auto cellValName : cellValNames)
    {
        this->c_d[cellValName].resize(msh.nCells);
    }
    for (auto planeValName : planeValNames)
    {
        this->p_d[planeValName].resize(msh.nPlanes);
    }
}
