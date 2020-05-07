//
// Created by Elijah on 03/03/2018.
//

#include "../../util/clUtils/clUtils.h"
#include "../../tests/run_tests/run_tests.h"
#include "../../util/simUtils/simUtils.h"
#include <CL/cl.h>
#include <malloc.h>
#include <string.h>

#ifndef DEMORANGES_SIMRUNNER_H
#define DEMORANGES_SIMRUNNER_H

int runSim(particle *hparticles, cl_ulong NUMPART, cl_kernel iterate_particle, cl_float particle_diameter,
           cl_bool periodic,
           cl_float domain_length, char prefix[], char log_dir[], float sim_length, float timestep,
           boolean VERBOSE, boolean LOG_DATA, boolean log_vel, float log_step,
           cl_device_id device, cl_context context, int coupled, int analytic, int fixed_Re);

#endif //DEMORANGES_SIMRUNNER_H
