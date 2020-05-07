//
// Created by Elijah on 05/12/2017.
//

#ifndef DEMORANGES_CLDEVICEUTILS_H
#define DEMORANGES_CLDEVICEUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <CL/cl.h>
#include "../../structures/particle.h"
#include <stdbool.h>
#include <malloc.h>
#if defined(_MSC_VER)
#include <Windows.h>
#endif

#define MAX_SOURCE_SIZE (0x100000)

// If boolean is not correctly defined.
#if !defined(boolean)
#define boolean bool
#if !defined(_MSC_VER)
#define TRUE true
#define FALSE false
#endif
#endif

void setContext(cl_device_id *device, cl_context *context,
                boolean verbose);

void printDeviceDetails(cl_uint platformCount, cl_platform_id *platforms);

cl_command_queue getCommandQueue(cl_context context, cl_device_id device, boolean verbose);

cl_context getContext(cl_device_id **devices, cl_uint num_devices, boolean verbose);

cl_kernel getKernelWithUtils(cl_device_id device, cl_context context, char *fileName, char *kernelName, boolean verbose);

cl_kernel
getKernel(cl_device_id device, cl_context context, char **fileNames, int numFiles, char *kernelName,
          boolean verbose);

cl_int particlesToDevice(cl_command_queue queue, cl_mem gparticles, particle **hparticles, cl_ulong NUMPART);

cl_int particlesToHost(cl_command_queue queue, cl_mem gparticles, particle **hparticles, cl_ulong NUMPART);

#endif //DEMORANGES_CLDEVICEUTILS_H
