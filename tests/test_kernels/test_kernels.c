//
// Created by Elijah on 19/12/2017.
//
#include "test_kernels.h"

boolean test_kernels(cl_device_id device, cl_context context, boolean verbose) {
    cl_int ret;

    if (verbose) printf("\nTesting compiling kernels.\n");

    // Source file variables
    char *fileNames[] = {
            PROJECT_DIR "/kernels/iterate_particle.cl"
    };

    char *kernelNames[] = {
            "iterate_particle"
    };
    int files = 1;

    for (int i = 0; i < files; i++) {

        if (verbose) {
            printf("\nKernel file: %s\n", fileNames[i]);
            printf("Kernel name: %s\n", kernelNames[i]);
        }

        cl_kernel kernel = getKernelWithUtils(device, context, fileNames[i], kernelNames[i], verbose);
        if (kernel == NULL) {
            return FALSE;
        }
    }
    return TRUE;
}
