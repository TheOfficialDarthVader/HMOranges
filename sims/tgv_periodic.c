#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else

#include <CL/cl.h>

#endif

#include <stdio.h>
#include <malloc.h>
#include "../util/clUtils/clUtils.h"
#include "../util/particleUtils/particleUtils.h"
#include "../tests/run_tests/run_tests.h"
#include "../util/simUtils/simUtils.h"
#include "simRunner/simRunner.h"

#define MAX_SOURCE_SIZE (0x100000)
#define VERBOSE TRUE
#define LOG_DATA FALSE

char *prefix = "TGV_PERIODIC";
char *log_dir = "TGV_PERIODIC/";

particle *hparticles;
cl_ulong NUMPART = 10000;

// Gas properties
cl_double C_p_G = 997;
cl_double W_G = 28.97; // Molecular weight of carrier gas (air) (kg/(kg mole))
cl_double Y_G = 0;
cl_double T_G = 298; // Gas temperature (K)
cl_double rho_G = 0.00001;

// Particle properties.
cl_double T_d = 282; //particle temperature (K)
cl_double density = 997;
cl_double particle_diameter = 0.0010488088481701515;
cl_double particle_effect_diameter;
cl_double fluid_viscosity;// = 1.04;//0.0000193
cl_double T_B = 373.15; // Boiling temperature
cl_double W_V = 18.015;  // Molecular weight of vapour phase (water) (kg/(kg mole))
cl_double C_L = 4184;  // (J/kg/K)

// Other properties
cl_double P_atm = 101325;  // Atmospheric pressure (Pa)
cl_double R_bar = 8314.5;
cl_double R = 287;  // Universal gas constant ()
cl_bool periodic = CL_TRUE;
cl_int coupled = 2;
cl_int analytic = 2;
cl_int fixed_Re = 0;

float init_speed_mean = 1;
float init_speed_std_dev = 0.1;

cl_float timestep;
cl_float sim_length;
cl_float log_step;//1.3847222222222226;

cl_float domain_length;

cl_context context;
cl_device_id device;

int main() {
    // Initialize OpenCL.
    setContext(&device, &context, TRUE);

    // Run tests
    if (!run_all_tests(device, context, FALSE)) {
        return 1;
    }

    // Build iterate_particle kernel.
    char *iterate_particle_files[] = {PROJECT_DIR "/util/kernelUtils.cl",
                                      PROJECT_DIR "/kernels/get_gravity/no_gravity.cl",
                                      PROJECT_DIR "/kernels/get_vel_fluid/tgv.cl",
                                      PROJECT_DIR "/kernels/iterate_particle.cl"};
    cl_kernel iterate_particle = getKernel(device, context, iterate_particle_files, 4, "iterate_particle", TRUE);

    hparticles = malloc(sizeof(particle) * NUMPART);
    if (hparticles == NULL) {
        fprintf(stderr, "Particles memory allocation failed.\n");
        return 1;
    }

    printf("[INIT] Creating particle positions.\n");
    particle_effect_diameter = (cl_float) (1.5 * particle_diameter);
    cl_float3 *positions = malloc(sizeof(cl_float3) * NUMPART);
    // Using particle_effect_diameter so that cohesion effects are considered at the appropriate range.
    float cube_length = createCubePositions(positions, NUMPART, particle_effect_diameter, 2, (cl_float3) {0, 0, 0});
    domain_length = (cl_float) (2 * PI);

    if (cube_length > domain_length) {
        fprintf(stderr, "Not all particles fit within the specified domain length for the given cube parameters (%.3f > %.3f).", cube_length, domain_length);
        return 1;
    }

    cl_float3 *velocities = malloc(sizeof(cl_float3) * NUMPART);
    createNormalDistVelocities(velocities, NUMPART, init_speed_mean, init_speed_std_dev);

    // Initialize particles.
    initializeMonodisperseParticles(hparticles, NUMPART, density, fluid_viscosity, particle_diameter,
                                    particle_effect_diameter, C_p_G, C_L, P_atm, rho_G, R_bar, R, W_G, W_V, Y_G, T_d, T_B, T_G, positions, velocities);
    free(positions);

    if (!checkPositions(hparticles, NUMPART, domain_length)) {
        fprintf(stderr, "Particles outside domain limits.\n");
        return 1;
    }

    float tau_d = (density*particle_diameter*particle_diameter)/(18*hparticles[0].fluid_viscosity);
    sim_length = 200.0*tau_d;
    timestep = (cl_float) (tau_d/8);

    if (!writeSetupData(prefix, log_dir, NUMPART, timestep, sim_length, domain_length,
                        particle_diameter, particle_effect_diameter, density, fluid_viscosity)) {
        return 1;
    }

    int sim_ret = runSim(hparticles, NUMPART, iterate_particle, particle_diameter, periodic, domain_length,
                         prefix, log_dir,
                         sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, device, context, coupled, analytic);

    clReleaseContext(context);
    return sim_ret;
}