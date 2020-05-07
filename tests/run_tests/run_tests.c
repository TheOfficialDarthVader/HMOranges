//
// Created by Elijah on 18/12/2017.
//

#include "run_tests.h"


boolean run_all_tests(cl_device_id device, cl_context context, boolean verbose){
    printf("Running tests.\n");

    if (!test_kernels(device, context, verbose)) {
        fprintf(stderr, "\nFAILED AT test_kernels\n");
        return FALSE;
    }

    if (!test_particle_struct_alignment(device, context, verbose)) {
        fprintf(stderr, "\nFAILED AT test_particle_struct_alignment\n");
        return FALSE;
    }
    printf("All tests passed.\n");

    return TRUE;
}

