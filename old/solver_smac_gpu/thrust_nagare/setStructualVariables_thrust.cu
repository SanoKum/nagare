#include "thrust_nagare/update_thrust.h"
#include "thrust_nagare/variables_thrust.h"

#include "thrust/device_vector.h"
#include "thrust/fill.h"

#include "flowFormat.hpp"

#include "variables.hpp"

void setStructualVariables_thrust(solverConfig& cfg , mesh& msh, variables_d& v_d)
{
    thrust::device_vector<flow_float>& Ux_d = v_d.c_d["Ux"];
    thrust::device_vector<flow_float>& Uy_d = v_d.c_d["Uy"];
    thrust::device_vector<flow_float>& Uz_d = v_d.c_d["Uz"];

    thrust::device_vector<flow_float>& T_d = v_d.c_d["T"];
    thrust::device_vector<flow_float>& P_d = v_d.c_d["T"];

    thrust::device_vector<flow_float>& divU_d = v_d.c_d["divU"];
    thrust::device_vector<flow_float>& cfl_d = v_d.c_d["cfl"];

}