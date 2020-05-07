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

#define MAX_SOURCE_SIZE (0x100000)
#define VERBOSE FALSE
#define LOG_DATA FALSE

char prefix[50];
char folder[100];
char dir[200];

particle *hparticles;
cl_ulong NUMPART = 1;

// Gas properties
cl_double C_p_G = 1007;
cl_double W_G = 28.97; // Molecular weight of carrier gas (air) (kg/(kg mole))
cl_double Y_G = 0;
cl_double T_G = 298; // Gas temperature (K)
cl_double rho_G = 1.184;

// Particle properties.
cl_double T_d = 282; //particle temperature (K)
cl_double density = 997;
cl_double particle_diameter;
cl_float particle_effect_diameter;
cl_double fluid_viscosity;// = 1.04;//0.0000193
cl_double T_B = 373.15; // Boiling temperature
cl_double W_V = 18.015;  // Molecular weight of vapour phase (water) (kg/(kg mole))
cl_double C_L = 4184;  // (J/kg/K)

// Other properties
cl_double P_atm = 101325;  // Atmospheric pressure (Pa)
cl_double R_bar = 8314.5;
cl_double R = 287;  // Universal gas constant ()
cl_bool periodic = CL_TRUE;
cl_int coupled = 1;
cl_int analytic = 2;
cl_int fixed_Re = 0;
double T_R;
float tau;
float Nu = 2;

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

    for (float i = 1; i <= 10; i += 1) {
        hparticles[0].diameter = powf(1.1,0.5)/1000;
        hparticles[0].effect_diameter = 0;
        hparticles[0].density = density;
        hparticles[0].id = 0;
        hparticles[0].pos = (cl_float3) {0, 0, 0};
        hparticles[0].vel = (cl_float3) {0, 0, 0};
        hparticles[0].forces = (cl_float3) {0, 0, 0};
        hparticles[0].fluid_vel = (cl_float3) {0.0, 0.0, 0.0};
        hparticles[0].m_d = density * M_PI * (particle_diameter * particle_diameter * particle_diameter) / 6;

        //Create reference temperature//
        T_R = (137 * powf((T_B/373.15), 0.68) * log10f(T_G)) - 45;

        //droplet/liquid properties
        hparticles[0].C_L = C_L;
        hparticles[0].T_d = T_d;
        hparticles[0].T_B = T_B;
        hparticles[0].Reynolds_d = 0;

        //vapour properties
        hparticles[0].L_V = (2.257*(pow(10,6)) + (2.595*(pow(10,3)) * (373.15-T_R)));
        hparticles[0].W_V = 18.015;

        //gas properties
        hparticles[0].P_atm = P_atm;
        hparticles[0].R_bar = R_bar;
        hparticles[0].R = R;
        hparticles[0].W_G = W_G;
        hparticles[0].Y_G = Y_G;
        hparticles[0].T_G = T_G;
        hparticles[0].rho_G = rho_G;
        hparticles[0].P_G = rho_G*R*T_R;
        hparticles[0].fluid_viscosity = (6.109*pow(10,-6)) + (4.604*pow(10,-8)*T_R) - (1.05*pow(10, -11)*T_R*T_R);

        //other properties
        hparticles[0].H_deltaT = 0;
        hparticles[0].f_2 = 1;
        hparticles[0].theta_1 = C_p_G/C_L;
        hparticles[0].theta_2 = W_G/W_V;

        //non-dimensional numbers
        hparticles[0].Pr_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R);
        hparticles[0].Sc_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R); // =Pr_G

        tau = get_tau(&(hparticles[0]));

        domain_length = 0.5;
        printf("Mass = %f\n", get_particle_mass(&(hparticles[0])));
        printf("Tau = %f\n", tau);

        timestep = tau * (i/10);
        log_step = timestep;
        sim_length = 6 * tau * ((3*hparticles[0].Pr_G)/(Nu*hparticles[0].f_2*hparticles[0].theta_1));
        printf("tau_h: %f\n", tau * ((3*hparticles[0].Pr_G)/(Nu*hparticles[0].f_2*hparticles[0].theta_1)));

        sprintf(prefix, "uc_heat_transfer");
        char c[50];
        float a = i/10.0;
        sprintf(c, "%.1f", a);
        char *token = strtok(c, ".");
        char *num[2];
        int j = 0;
        while (token){
            num[j++] = token;
            token = strtok(NULL, ".");
        }

        snprintf(folder, sizeof(folder), "%s%s%s%s%s", "uc_heat_transfer_", num[0], "_", num[1], "_tau/");
        sprintf(dir, "%s%s", PROJECT_DIR "/verification/uncoupled_heat_transfer/data/", folder);

        runSim(hparticles, NUMPART, iterate_particle, particle_diameter, periodic, domain_length,
               prefix, dir,
               sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, device, context, coupled, analytic, fixed_Re);
    }

//    timestep = tau * 0.001;
//    log_step = timestep;
//
//    snprintf(folder, sizeof(folder), "%s%s%s%s%s", "uc_heat_transfer_", "0", "_", "001", "_tau/");
//    sprintf(dir, "%s%s", PROJECT_DIR "/verification/uncoupled_heat_transfer/data/", folder);
//    hparticles[0].diameter = powf(1.1,power_half)/1000;
//    hparticles[0].effect_diameter = 0;
//    hparticles[0].fluid_viscosity = fluid_viscosity;
//    hparticles[0].density = density;
//    hparticles[0].id = 0;
//    hparticles[0].pos = (cl_float3) {0, 0, 0};
//    hparticles[0].vel = (cl_float3) {0, 0, 0};
//    hparticles[0].forces = (cl_float3) {0, 0, 0};
//    hparticles[0].fluid_vel = (cl_float3) {0.0, 0.0, 0.0};
//    hparticles[0].m_d = density * M_PI * (hparticles[0].diameter * hparticles[0].diameter * hparticles[0].diameter) / 6;
//
//    //Create reference temperature//
//    T_R = (137 * powf((T_B/373.15), 0.68) * log10f(T_G)) - 45;
//
//    //droplet/liquid properties
//    hparticles[0].C_L = C_L;
//    hparticles[0].T_d = T_d;
//    hparticles[0].T_B = T_B;
//
//    //vapour properties
//    hparticles[0].L_V = (2.257*(pow(10,6)) + (2.595*(pow(10,3)) * (373.15-T_R)));
//    hparticles[0].W_V = 18.015;
//
//    //gas properties
//    hparticles[0].P_atm = P_atm;
//    hparticles[0].R_bar = R_bar;
//    hparticles[0].R = R;
//    hparticles[0].W_G = W_G;
//    hparticles[0].Y_G = Y_G;
//    hparticles[0].T_G = T_G;
//    hparticles[0].P_G = rho_G*R*T_R;
//    hparticles[0].fluid_viscosity = (6.109*pow(10,-6)) + (4.604*pow(10,-8)*T_R) - (1.05*pow(10, -11)*T_R*T_R);
//
//    //other properties
//    hparticles[0].H_deltaT = 0;
//    hparticles[0].f_2 = 1;
//    hparticles[0].theta_1 = C_p_G/C_L;
//    hparticles[0].theta_2 = W_G/W_V;
//
//    //non-dimensional numbers
//    hparticles[0].Pr_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R);
//    hparticles[0].Sc_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R); // =Pr_G
//
//    runSim(hparticles, NUMPART, iterate_particle, particle_diameter, periodic, domain_length,
//           prefix, dir,
//           sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, device, context, coupled, analytic);
//
//    //analytic solution
//    timestep = tau * 0.1;
//    log_step = timestep;
//    analytic = 1;
//
//    sprintf(prefix, "uc_an_heat_transfer");
//    snprintf(folder, sizeof(folder), "%s%s%s%s%s", "an_uc_heat_transfer_", "0", "_", "1", "_tau/");
//    sprintf(dir, "%s%s", PROJECT_DIR "/verification/uncoupled_heat_transfer/data/", folder);
//
//    hparticles[0].diameter = powf(1.1,power_half)/1000;
//    hparticles[0].effect_diameter = 0;
//    hparticles[0].fluid_viscosity = fluid_viscosity;
//    hparticles[0].density = density;
//    hparticles[0].id = 0;
//    hparticles[0].pos = (cl_float3) {0, 0, 0};
//    hparticles[0].vel = (cl_float3) {0, 0, 0};
//    hparticles[0].forces = (cl_float3) {0, 0, 0};
//    hparticles[0].fluid_vel = (cl_float3) {0.0, 0.0, 0.0};
//    hparticles[0].m_d = density * M_PI * (particle_diameter * particle_diameter * particle_diameter) / 6;
//
//    //Create reference temperature//
//    T_R = (137 * powf((T_B/373.15), 0.68) * log10f(T_G)) - 45;
//
//    //droplet/liquid properties
//    hparticles[0].C_L = C_L;
//    hparticles[0].T_d = T_d;
//    hparticles[0].T_B = T_B;
//
//    //vapour properties
//    hparticles[0].L_V = (2.257*(pow(10,6)) + (2.595*(pow(10,3)) * (373.15-T_R)));
//    hparticles[0].W_V = 18.015;
//
//    //gas properties
//    hparticles[0].P_atm = P_atm;
//    hparticles[0].R_bar = R_bar;
//    hparticles[0].R = R;
//    hparticles[0].W_G = W_G;
//    hparticles[0].Y_G = Y_G;
//    hparticles[0].T_G = T_G;
//    hparticles[0].P_G = rho_G*R*T_R;
//    hparticles[0].fluid_viscosity = (6.109*pow(10,-6)) + (4.604*pow(10,-8)*T_R) - (1.05*pow(10, -11)*T_R*T_R);
//
//    //other properties
//    hparticles[0].H_deltaT = 0;
//    hparticles[0].f_2 = 1;
//    hparticles[0].theta_1 = C_p_G/C_L;
//    hparticles[0].theta_2 = W_G/W_V;
//
//    //non-dimensional numbers
//    hparticles[0].Pr_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R);
//    hparticles[0].Sc_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R); // =Pr_G
//
//    runSim(hparticles, NUMPART, iterate_particle, particle_diameter, periodic, domain_length,
//           prefix, dir,
//           sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, device, context, coupled, analytic, fixed_Re);

}