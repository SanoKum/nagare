#include "thrust_nagare/update_thrust.h"
#include "thrust_nagare/variables_thrust.h"

#include "thrust/device_vector.h"
#include "thrust/fill.h"

#include "flowFormat.hpp"

void updateVariablesForNextLoop_thrust(solverConfig& cfg, mesh& msh, variables_d& v_d, matrix& mat_p)
{
    thrust::device_vector<flow_float>& US = v_d.p_d["US"];
    thrust::device_vector<flow_float>& USN = v_d.p_d["USN"];

    thrust::device_vector<flow_float>& Ux = v_d.c_d["Ux"];
    thrust::device_vector<flow_float>& Uy = v_d.c_d["Uy"];
    thrust::device_vector<flow_float>& Uz = v_d.c_d["Uz"];

    thrust::device_vector<flow_float>& UxN = v_d.c_d["UxN"];
    thrust::device_vector<flow_float>& UyN = v_d.c_d["UyN"];
    thrust::device_vector<flow_float>& UzN = v_d.c_d["UzN"];

    thrust::device_vector<flow_float>& T = v_d.c_d["T"];
    thrust::device_vector<flow_float>& TN = v_d.c_d["TN"];

    thrust::device_vector<flow_float>& dUxdx = v_d.c_d["dUxdx"];
    thrust::device_vector<flow_float>& dUydx = v_d.c_d["dUydx"];
    thrust::device_vector<flow_float>& dUzdx = v_d.c_d["dUzdx"];

    thrust::device_vector<flow_float>& dUxdy = v_d.c_d["dUxdy"];
    thrust::device_vector<flow_float>& dUydy = v_d.c_d["dUydy"];
    thrust::device_vector<flow_float>& dUzdy = v_d.c_d["dUzdy"];

    thrust::device_vector<flow_float>& dUxdz = v_d.c_d["dUxdz"];
    thrust::device_vector<flow_float>& dUydz = v_d.c_d["dUydz"];
    thrust::device_vector<flow_float>& dUzdz = v_d.c_d["dUzdz"];

    thrust::device_vector<flow_float>& dPdx = v_d.c_d["dPdx"];
    thrust::device_vector<flow_float>& dPdy = v_d.c_d["dPdy"];
    thrust::device_vector<flow_float>& dPdz = v_d.c_d["dPdz"];

    thrust::device_vector<flow_float>& convx = v_d.c_d["convx"];
    thrust::device_vector<flow_float>& convy = v_d.c_d["convy"];
    thrust::device_vector<flow_float>& convz = v_d.c_d["convz"];

    thrust::device_vector<flow_float>& diffx = v_d.c_d["diffx"];
    thrust::device_vector<flow_float>& diffy = v_d.c_d["diffy"];
    thrust::device_vector<flow_float>& diffz = v_d.c_d["diffz"];

    thrust::device_vector<flow_float>& src = v_d.c_d["src"];
    thrust::device_vector<flow_float>& res = v_d.c_d["res"];

    thrust::device_vector<flow_float>& diffT = v_d.c_d["diffT"];
    thrust::device_vector<flow_float>& dP = v_d.c_d["dP"];
    thrust::device_vector<flow_float>& dPPdx = v_d.c_d["dPPdx"];
    thrust::device_vector<flow_float>& dPPdy = v_d.c_d["dPPdy"];
    thrust::device_vector<flow_float>& dPPdz = v_d.c_d["dPPdz"];
    thrust::device_vector<flow_float>& divU = v_d.c_d["divU"];
    thrust::device_vector<flow_float>& divU_vol = v_d.c_d["divU*vol"];
    thrust::device_vector<flow_float>& divU_star = v_d.c_d["divU_star"];

    // d -> d
    thrust::copy(Ux.begin(), Ux.end(), UxN.begin());
    thrust::copy(Uy.begin(), Uy.end(), UyN.begin());
    thrust::copy(Uz.begin(), Uz.end(), UzN.begin());
    thrust::copy(T.begin(), T.end(), TN.begin());


    thrust::fill(convx.begin(), convx.end(), 0.0);
    thrust::fill(convy.begin(), convy.end(), 0.0);
    thrust::fill(convz.begin(), convz.end(), 0.0);

    thrust::fill(diffx.begin(), diffx.end(), 0.0);
    thrust::fill(diffy.begin(), diffy.end(), 0.0);
    thrust::fill(diffz.begin(), diffz.end(), 0.0);

    thrust::fill(dP.begin(), dP.end(), 0.0);

    thrust::fill(dPPdx.begin(), dPPdx.end(), 0.0);
    thrust::fill(dPPdy.begin(), dPPdy.end(), 0.0);
    thrust::fill(dPPdz.begin(), dPPdz.end(), 0.0);

    thrust::fill(divU.begin(), divU.end(), 0.0);
    thrust::fill(divU_vol.begin(), divU_vol.end(), 0.0);
    thrust::fill(divU_star.begin(), divU_star.end(), 0.0);

}