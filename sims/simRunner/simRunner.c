//
// Created by Elijah on 03/03/2018.
//

#include "simRunner.h"


int runSim(particle *hparticles, cl_ulong NUMPART, cl_kernel iterate_particle, cl_float particle_diameter,
           cl_bool periodic,
           cl_float domain_length, char prefix[], char log_dir[], float sim_length, float timestep,
           boolean VERBOSE, boolean LOG_DATA, boolean log_vel, float log_step,
           cl_device_id device, cl_context context, cl_int coupled, cl_int analytic, cl_int fixed_Re) {
    cl_int ret;  // Variable in which to store OpenCL return values.

    cl_mem gparticles;  // GPU array of particles.
    
    // Check that the logging directory exists if needed.
    if (!checkDirExists(log_dir) && LOG_DATA) {
        fprintf(stderr, "Error: Directory (%s) does not exist or cannot be accessed.\n", log_dir);
        return 1;
    }

    // Create command queue.
    cl_command_queue queue = getCommandQueue(context, device, VERBOSE);

    // Create gparticles buffer and copy particles into it.
    gparticles = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(particle) * NUMPART, NULL, &ret);
    ret = particlesToDevice(queue, gparticles, &hparticles, NUMPART);

    int int_periodic;  // Enumerated version of periodic boolean as booleans can't be passed to GPU by default.
    if (periodic) {
        int_periodic = 1;
    } else {
        int_periodic = 0;
    }

    // Set all kernel arguments. Note that if the cl_mem object changes they must be reset (e.g. gpp_cols).
    printf("[INIT] Setting kernel arguments\n");
    ret = clSetKernelArg(iterate_particle, 0, sizeof(cl_mem), &gparticles);
    ret = clSetKernelArg(iterate_particle, 1, sizeof(cl_float), &timestep);
    ret = clSetKernelArg(iterate_particle, 2, sizeof(cl_int), &int_periodic);
    ret = clSetKernelArg(iterate_particle, 3, sizeof(cl_float), &domain_length);
    ret = clSetKernelArg(iterate_particle, 4, sizeof(cl_int), &coupled);
    ret = clSetKernelArg(iterate_particle, 5, sizeof(cl_int), &analytic);
    ret = clSetKernelArg(iterate_particle, 6, sizeof(cl_int), &fixed_Re);

    // Do pre-simulation output and logging.
    printf("\nRunning sim with %llu particles, timestep %f, and log step %f, length %f.\n", NUMPART, timestep, log_step,
           sim_length);

    if (!writeTime(prefix, log_dir, NUMPART, "Start")) {
        return 1;
    }
    if (LOG_DATA) {
        printf("Logging at time: 0.000000\n");
        if (!writeParticles(hparticles, 0, prefix, log_dir, NUMPART, log_vel)) {
            return 1;
        }
    }

    cl_float last_write = 0; // Variable to store when logging was last run.

    log_step = timestep;

    // Run simulation.
    for (cl_float time = timestep; time <= sim_length; time += timestep) {
        if (VERBOSE) printf("Time = %f\n", time);

        if (ret != 0) {
            fprintf(stderr, "Error: Failed to read from collision count buffer. Error code %i\n", ret);
            return 1;
        }
        // Iterate particles.
        if (VERBOSE) printf("    Iterating particles\n");
        ret = clEnqueueNDRangeKernel(queue, iterate_particle, 1, NULL, &NUMPART, 0, NULL, NULL, NULL);

        if (LOG_DATA && time - last_write >= 0.99*log_step) {
            ret = particlesToHost(queue, gparticles, &hparticles, NUMPART);
            printf("Logging at time: %f\n", time);

            if (!writeParticles(hparticles, time, prefix, log_dir, NUMPART, log_vel)) {
                return 1;
            }

            last_write = time;
        }
    }

//    if (LOG_DATA) {
//        printf("Logging at time: %f\n", sim_length);
//        if (!writeParticles(hparticles, sim_length, prefix, log_dir, NUMPART, log_vel)) {
//            return 1;
//        }
//    }

    if (!writeTime(prefix, log_dir, NUMPART, "End")) {
        return 1;
    }
    return 0;
}