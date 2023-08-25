#include "thrust_nagare/setInitial_thrust.h"

#include "thrust/device_vector.h"
#include "thrust/fill.h"

#include "flowFormat.hpp"

void setInitial_thrust(solverConfig& cfg , mesh& msh , variables_d& v)
{
    thrust::device_vector<flow_float>& T = v.c_d["T"];
    thrust::device_vector<flow_float>& Ux= v.c_d["Ux"];
    thrust::device_vector<flow_float>& Uy= v.c_d["Uy"];
    thrust::device_vector<flow_float>& Uz= v.c_d["Uy"];
    thrust::device_vector<flow_float>& P = v.c_d["P"];
    thrust::device_vector<flow_float>& ro= v.c_d["ro"];

    thrust::fill(T.begin() , T.end(), 300.0);
    thrust::fill(Ux.begin(), Ux.end(), 5.0);
    thrust::fill(Uy.begin(), Uy.end(), 5.0);
    thrust::fill(Uz.begin(), Uz.end(), 5.0);
    thrust::fill(P.begin() , P.end() , 5.0);
    thrust::fill(ro.begin(), ro.end(), 1.2);

    std::cout << "filled" << std::endl;
};