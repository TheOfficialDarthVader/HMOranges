//
// Created by Elijah on 19/12/2017.
//

#ifndef DEMORANGES_TEST_ALIGNMENT_H
#define DEMORANGES_TEST_ALIGNMENT_H

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <CL/cl.h>
#include "../../util/clUtils/clUtils.h"
#include "../../structures/particle.h"

boolean test_particle_struct_alignment(cl_device_id device, cl_context context, boolean verbose);

#endif //DEMORANGES_TEST_ALIGNMENT_H
