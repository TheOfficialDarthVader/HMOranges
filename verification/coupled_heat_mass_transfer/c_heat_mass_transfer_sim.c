//
// Created by akern on 02/04/2020.
//

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else

#include <CL/cl.h>

#endif

#include <stdio.h>
#include "../../util/clUtils/clUtils.h"
#include "../../util/particleUtils/particleUtils.h"
#include "../../tests/run_tests/run_tests.h"
#include "../../util/simUtils/simUtils.h"
#include "../../sims/simRunner/simRunner.h"
#include <malloc.h>
#include <string.h>

#define MAX_SOURCE_SIZE (0x100000)
#define VERBOSE FALSE
#define LOG_DATA FALSE

char prefix[50];
char folder[100];
char dir[200];

particle *hparticles;
cl_ulong NUMPART = 1;

// Gas properties
cl_double fluid_viscosity;// = 1.04;//0.0000193
cl_double W_G = 28.97; // Molecular weight of carrier gas (air) (kg/(kg mole))
cl_double Y_G = 0.0;

// Particle properties.
cl_float particle_effect_diameter;

// Other properties
cl_double P_atm = 101325.0;  // Atmospheric pressure (Pa)
cl_double R_bar = 8314.5;
cl_double R = 287.0;  // Universal gas constant ()
cl_bool periodic = CL_TRUE;
cl_int coupled = 2;
cl_int analytic = 2;
cl_int fixed_Re = 1;
double T_R;
float tau;

cl_float timestep;
cl_float sim_length;
cl_float log_step;

cl_float domain_length;

cl_context context;
cl_device_id device;

int main() {

    // Initializing OpenCL.
    setContext(&device, &context, TRUE);

    // Run tests
    if (!run_all_tests(device, context, FALSE)) {
        return 1;
    }

    char *iterate_particle_files[] = {PROJECT_DIR "/util/kernelUtils.cl",
                                      PROJECT_DIR "/kernels/get_gravity/gravity.cl",
                                      PROJECT_DIR "/kernels/get_vel_fluid/no_fluid_flow.cl",
                                      PROJECT_DIR "/kernels/iterate_particle.cl"};
    cl_kernel iterate_particle = getKernel(device, context, iterate_particle_files, 4, "iterate_particle", TRUE);

    hparticles = malloc(sizeof(particle) * NUMPART);

    printf("%u\n", sizeof(P_atm));
    printf("%u\n", sizeof(particle));

    char drop_species[3][50] = {"water", "hexane", "decane"};

    float sim_length_vals[3] = {260.0, 2.0, 1.0};

    double density_vals[3] = {997.0, 664.0, 642.0};
    double T_B_vals[3] = {373.15, 344.6, 447.7};
    double W_V_vals[3] = {18.015, 86.178, 142.0};
    double C_L_vals[3] = {4184.0, 2302.0, 2520.5};

    double T_G_vals[3] = {298.0, 437.0, 1000.0};
    double rho_G_vals[3] = {1.184, 0.807, 0.3529};
    double C_p_G_vals[3] = {1007.0, 1020.0, 1141.0};
    double Re_d_vals[3] = {0.0, 110.0, 17.0};
    double T_d_vals[3] = {282.0, 281.0, 315.0};
    double D_vals[3] = {pow(1.1,0.5)/1000, 0.00176, 0.002};

//    for (int i = 0; i <= 2; i += 1) {
//        printf("%s\n", drop_species[i]);
//
//        hparticles[0].diameter = D_vals[i];
//        hparticles[0].effect_diameter = 0;
//        hparticles[0].density = density_vals[i];
//        hparticles[0].id = 0;
//        hparticles[0].pos = (cl_float3) {0, 0, 0};
//        hparticles[0].vel = (cl_float3) {0, 0, 0};
//        hparticles[0].forces = (cl_float3) {0, 0, 0};
//        hparticles[0].fluid_vel = (cl_float3) {0.0, 0.0, 0.0};
//        hparticles[0].m_d = hparticles[0].density * M_PI * (hparticles[0].diameter * hparticles[0].diameter * hparticles[0].diameter) / 6;
//        hparticles[0].initial_mass = hparticles[0].density * M_PI * (hparticles[0].diameter * hparticles[0].diameter * hparticles[0].diameter) / 6;
//
//        //Create reference temperature//
//        T_R = (137 * pow((T_B_vals[i]/373.15), 0.68) * log10(T_G_vals[i])) - 45;
//
//        //droplet/liquid properties
//        hparticles[0].C_L = C_L_vals[i];
//        hparticles[0].T_d = T_d_vals[i];
//        hparticles[0].T_B = T_B_vals[i];
//        hparticles[0].Reynolds_d = Re_d_vals[i];
//
//        //vapour properties
//        if (strcmp(drop_species[i], "water") == 0) {
//            hparticles[0].L_V = (2.257*(pow(10,6)) + (2.595*(pow(10,3)) * (373.15-T_R)));
//        } if (strcmp(drop_species[i], "hexane") == 0) {
//            hparticles[0].L_V = (5.1478*pow(10,5))*pow((1-(T_R/512.0)),0.3861);
//        } if (strcmp(drop_species[i], "decane") == 0) {
//            hparticles[0].L_V = (3.958*pow(10,4))*(pow(619.0-T_R, 0.38));
//            printf("%.10lf\n", (3.958*pow(10,4))*(pow(619.0-T_R, 0.38)));
//        }
//
//        hparticles[0].W_V = W_V_vals[i];
//
//        //gas properties
//        hparticles[0].P_atm = P_atm;
//        hparticles[0].R_bar = R_bar;
//        hparticles[0].R = R;
//        hparticles[0].W_G = W_G;
//        hparticles[0].Y_G = Y_G;
//        hparticles[0].T_G = T_G_vals[i];
//        hparticles[0].rho_G = rho_G_vals[i];
//        hparticles[0].P_G = rho_G_vals[i]*R*T_R;
//        hparticles[0].fluid_viscosity = (6.109*pow(10,-6)) + (4.604*pow(10,-8)*T_R) - (1.05*pow(10, -11)*T_R*T_R);
//
//        //other properties
//        hparticles[0].H_deltaT = 0;
//        hparticles[0].f_2 = 1;
//        hparticles[0].theta_1 = C_p_G_vals[i]/C_L_vals[i];
//        hparticles[0].theta_2 = W_G/W_V_vals[i];
//
//        //non-dimensional numbers
//        hparticles[0].Pr_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R);
//        hparticles[0].Sc_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R); // =Pr_G
//
//        tau = get_tau(&(hparticles[0]));
//
//        domain_length = 0.5;
//        printf("Mass = %f\n", get_particle_mass(&(hparticles[0])));
//        printf("Tau = %f\n", tau);
//
//        timestep = tau / (64);
//        log_step = timestep;
//        sim_length =  sim_length_vals[i] * tau;
//
//        snprintf(folder, sizeof(folder), "%s%s%s", "c_heat_mass_transfer_", drop_species[i], "/");
//        printf("%s\n", folder);
//        sprintf(dir, "%s%s", PROJECT_DIR "/verification/coupled_heat_mass_transfer/data/", folder);
//
//        runSim(hparticles, NUMPART, iterate_particle, hparticles[0].diameter, periodic, domain_length,
//               prefix, dir,
//               sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, device, context, coupled, analytic, fixed_Re);
//    }

}