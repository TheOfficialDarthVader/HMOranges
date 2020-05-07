//
// Created by Elijah on 30/11/2017.
//

#ifndef DEMORANGES_PARTICLE_H
#define DEMORANGES_PARTICLE_H

#include <CL/cl_platform.h>

typedef struct particle {
    cl_float3 pos;
    cl_float3 vel;
    cl_float3 forces;
    cl_float3 fluid_vel;
    cl_ulong id;
    cl_ulong cv_array_idx;
    cl_double diameter;
    cl_float effect_diameter;
    cl_double density;
    cl_double fluid_viscosity;
    cl_double T_d;
    cl_double W_V;
    cl_double T_B;
    cl_double L_V;
    cl_double C_L;
    cl_double P_atm;
    cl_double R_bar;
    cl_double R;
    cl_double W_G;
    cl_double theta_1;
    cl_double theta_2;
    cl_double Y_G;
    cl_double Pr_G;
    cl_double Sc_G;
    cl_double f_2;
    cl_double P_G;
    cl_double T_G;
    cl_double rho_G;
    cl_double H_deltaT;
    cl_double m_d;
    cl_double Reynolds_d;
    cl_double initial_mass;
    // Structure memory alignment for Visual Studio and GCC compilers.
#if defined(_MSC_VER)
    cl_char padding[244];
} __declspec(align(512)) particle;
#elif defined(__GNUC__) || defined(__GNUG__) || defined(__MINGW_GCC_VERSION)
} __attribute__((aligned (512))) particle;
#endif

#endif //DEMORANGES_PARTICLE_H
