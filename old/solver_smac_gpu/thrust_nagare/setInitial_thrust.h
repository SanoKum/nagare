#pragma once

#include "mesh/mesh.hpp"
#include "thrust_nagare/variables_thrust.h"
#include "input/solverConfig.hpp"

#include <cmath>

void setInitial_thrust(solverConfig& , mesh& , variables_d& );
