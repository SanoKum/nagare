#pragma once

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"

#include "variables.hpp"
#include "thrust_nagare/variables_thrust.h"

//void outputH5_XDMF(const mesh &msh , const variables &var)
void setStructualVariables_thrust(solverConfig& , mesh& , variables_d&);