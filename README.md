# HMOranges
HMOranges builds off [Oranges](https://github.com/TheOfficialDarthVader/Oranges-master) by adding a heat and mass transfer model for the droplets. The implementation is as per Oranges, with extra variables and a change to the kernel to include the heat and mass transfer.

## Project Structure
The project is organised into relevant directories. The contents of each of these is described below.

#### analysis/
This directory contains all of the analysis tools used to post-process data from simulations.
The primary functionality is the `Final_Stats_Test.py` which is used to collate and graph the data from simulations.

#### kernels/
This directory contains all of the OpenCL kernels used in the project.
When running the project this directory must be in the same directory as the working directory so that the program can access the kernel files.
Note, the kernel utilities file is in `util`, not `kernels`.

####  sims/
This directory contains all simulation programs as well as the main `simRunner` files.
`simRunner` is accessed by the simulations and handles all of the simulation running code.
A variety of simulations are available and are all CMake targets.

#### structures/
This directory contains all of the data structures used in the host code.
Matching device structures are found in `util/kernelUtils.cl`.

#### tests/
This directory contains all of the unit tests for the project.
These are accessed from `test/run_tests/` in simulations or can be run separately with more debug messages with the `run_tests` target.

#### util/
This directory contains all of the utility functions used in the project.
Most of these are found in `util/clUtils/` which contains all of the OpenCL utility functions that make calling OpenCL functions much simpler.
This directory also contains `kernelUtils.cl` which contains all of the utility functions and data structures for the OpenCL kernels.

#### verification/
This directory contains all of the verification test cases.
Actual simulations can be run using their targets and the resulting data can be graphed against the analytic solutions with the Python script in each directory.
