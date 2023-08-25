#include "thrust_nagare/update_thrust.h"
#include "thrust_nagare/variables_thrust.h"

#include "thrust/device_vector.h"
#include "thrust/fill.h"

#include "flowFormat.hpp"

#include "variables.hpp"

void copyDeviceToHost_cell(variables_d& v_d , variables& v)
{
    thrust::device_vector<flow_float>& Ux_d = v_d.c_d["Ux"];
    thrust::device_vector<flow_float>& Uy_d = v_d.c_d["Uy"];
    thrust::device_vector<flow_float>& Uz_d = v_d.c_d["Uz"];

    thrust::device_vector<flow_float>& T_d = v_d.c_d["T"];
    thrust::device_vector<flow_float>& P_d = v_d.c_d["T"];

    thrust::device_vector<flow_float>& divU_d = v_d.c_d["divU"];
    thrust::device_vector<flow_float>& cfl_d = v_d.c_d["cfl"];


    //thrust::host_vector<flow_float>& Ux = v.c["Ux"];
    //thrust::host_vector<flow_float>& Uy = v.c["Uy"];
    //thrust::host_vector<flow_float>& Uz = v.c["Uz"];

    //thrust::host_vector<flow_float>& T = v.c["T"];
    //thrust::host_vector<flow_float>& P = v.c["T"];

    //thrust::host_vector<flow_float>& divU = v.c["divU"];


    std::vector<flow_float>& Ux = v.c["Ux"];
    std::vector<flow_float>& Uy = v.c["Uy"];
    std::vector<flow_float>& Uz = v.c["Uz"];

    std::vector<flow_float>& T = v.c["T"];
    std::vector<flow_float>& P = v.c["P"];

    std::vector<flow_float>& divU = v.c["divU"];
    std::vector<flow_float>& cfl = v.c["cfl"];

    thrust::copy(Ux_d.begin(), Ux_d.end(), Ux.begin());
    thrust::copy(Uy_d.begin(), Uy_d.end(), Uy.begin());
    thrust::copy(Uz_d.begin(), Uz_d.end(), Uz.begin());

    thrust::copy(T_d.begin(), T_d.end(), T.begin());
    thrust::copy(P_d.begin(), P_d.end(), P.begin());
    thrust::copy(divU_d.begin(), divU_d.end(), divU.begin());
    thrust::copy(cfl_d.begin(), cfl_d.end(), cfl.begin());

}