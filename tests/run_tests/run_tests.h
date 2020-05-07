//
// Created by Elijah on 14/02/2018.
//

#ifndef DEMORANGES_RUN_TESTS_H
#define DEMORANGES_RUN_TESTS_H

#include "../test_kernels/test_kernels.h"
#include "../test_alignment/test_alignment.h"
#include "../test_setContext/test_setContext.h"

boolean run_all_tests(cl_device_id device, cl_context context, boolean verbose);

#endif //DEMORANGES_RUN_TESTS_H
