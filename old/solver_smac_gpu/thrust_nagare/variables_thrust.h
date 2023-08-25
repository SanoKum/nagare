#pragma once

#include <iostream>

#include<map>
#include"flowFormat.hpp"

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/sort.h>

#include <list>
#include <vector>
#include <string>

class variables_d {
public:
    std::map<std::string, thrust::device_vector<flow_float>> c_d;
    std::map<std::string, thrust::device_vector<flow_float>> p_d;

    thrust::device_vector<thrust::device_vector<geom_int>> pcmap;

    variables_d(mesh& , std::list<std::string> , std::list<std::string>);
};