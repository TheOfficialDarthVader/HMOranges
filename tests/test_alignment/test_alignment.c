//
// Created by Elijah on 19/12/2017.
//

#include "test_alignment.h"

boolean test_particle_struct_alignment(cl_device_id device, cl_context context, boolean verbose) {
    particle *hparticles;
    cl_mem gparticles;
    cl_ulong NUMPART = 10;

    cl_int ret;
    cl_bool hcorrect = TRUE;
    cl_mem gcorrect;

    if (verbose) printf("\nTesting particle struct alignment.\n");

    cl_kernel kernel = getKernelWithUtils(device, context, PROJECT_DIR "/tests/test_alignment/alignment_test_kernels.cl",
                                          "test_particle_struct_alignment", verbose);
    cl_command_queue queue = getCommandQueue(context, device, verbose);

    hparticles = malloc(sizeof(particle) * NUMPART);

    for (cl_ulong i = 0; i < NUMPART; i++) {
        hparticles[i].pos = (cl_float3) {20, 21, 22};
        hparticles[i].vel = (cl_float3) {23, 24, 25};
        hparticles[i].forces = (cl_float3) {26, 27, 28};
        hparticles[i].fluid_vel = (cl_float3) {29, 30, 31};
        hparticles[i].id = 32;
        hparticles[i].cv_array_idx = 33;
        hparticles[i].diameter = 34;
        hparticles[i].effect_diameter = 35;
        hparticles[i].density = 36;
        hparticles[i].fluid_viscosity = 37;
        hparticles[i].T_d = 38;
        hparticles[i].W_V = 39;
        hparticles[i].T_B = 40;
        hparticles[i].L_V = 41;
        hparticles[i].C_L = 42;
        hparticles[i].P_atm = 43;
        hparticles[i].R_bar = 44;
        hparticles[i].R = 45;
        hparticles[i].W_G = 46;
        hparticles[i].theta_1 = 47;
        hparticles[i].theta_2 = 48;
        hparticles[i].Y_G = 49;
        hparticles[i].Pr_G = 50;
        hparticles[i].Sc_G = 51;
        hparticles[i].f_2 = 52;
        hparticles[i].P_G = 53;
        hparticles[i].T_G = 54;
        hparticles[i].H_deltaT = 55;
        hparticles[i].m_d = 56;
    }

    gparticles = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(particle) * NUMPART, NULL, &ret);
    gcorrect = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(boolean), NULL, &ret);

    ret = particlesToDevice(queue, gparticles, &hparticles, NUMPART);
    ret = clEnqueueWriteBuffer(queue, gcorrect, CL_TRUE, 0, sizeof(boolean), &hcorrect, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), &gparticles);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), &gcorrect);

    ret = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, (size_t *) &NUMPART, 0, NULL, NULL, NULL);

    ret = particlesToHost(queue, gparticles, &hparticles, NUMPART);
    ret = clEnqueueReadBuffer(queue, gcorrect, CL_TRUE, 0, sizeof(boolean), &hcorrect, 0, NULL, NULL);

    ret = clFinish(queue);

    if (!hcorrect) {
        return FALSE;
    }

    for (int i = 0; i < NUMPART; i++) {
        particle p = hparticles[i];
        if (!(p.pos.x == 60 &&
              p.pos.y == 61 &&
              p.pos.z == 62 &&
              p.vel.x == 63 &&
              p.vel.y == 64 &&
              p.vel.z == 65 &&
              p.forces.x == 66 &&
              p.forces.y == 67 &&
              p.forces.z == 68 &&
              p.fluid_vel.x == 69 &&
              p.fluid_vel.y == 70 &&
              p.fluid_vel.z == 71 &&
              p.id == 72 &&
              p.cv_array_idx == 73 &&
              p.diameter == 74 &&
              p.effect_diameter == 75 &&
              p.density == 76 &&
              p.fluid_viscosity == 77 &&
              p.T_d == 78 &&
              p.W_V == 79 &&
              p.T_B == 80 &&
              p.L_V == 81 &&
              p.C_L == 82 &&
              p.P_atm == 83 &&
              p.R_bar == 84 &&
              p.R == 85 &&
              p.W_G == 86 &&
              p.theta_1 == 87 &&
              p.theta_2 == 88 &&
              p.Y_G == 89 &&
              p.Pr_G == 90 &&
              p.Sc_G == 91 &&
              p.f_2 == 92 &&
              p.P_G == 93 &&
              p.T_G == 94 &&
              p.H_deltaT == 95 &&
              p.m_d == 96
              )) {
            hcorrect = FALSE;
        }
    }
    return (boolean) hcorrect;
}
