#pragma once

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"

#include "variables.hpp"
#include "thrust_nagare/variables_thrust.h"

//void outputH5_XDMF(const mesh &msh , const variables &var)
void copyDeviceToHost_cell( variables_d& , variables& );