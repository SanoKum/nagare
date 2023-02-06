#pragma once
#include <iostream>
#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"

#include "variables.hpp"

//TODO: use template???

namespace cudaWrapper{
    void cudaMalloc_wrapper(flow_float* , geom_int );
    void cudaMalloc_wrapper(geom_int* , geom_int );

    void cudaMemcpy_vectorToDevice_wrapper(std::vector<flow_float>& , flow_float* );
    void cudaMemcpy_deviceToVector_wrapper(flow_float* , std::vector<flow_float>& );

    void cudaFree_wrapper(flow_float* );
    void cudaFree_wrapper(geom_int* );

    //void copyVariables_cell_plane_H2D(variables& );
}
