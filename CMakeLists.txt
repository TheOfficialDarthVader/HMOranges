cmake_minimum_required(VERSION 3.7)
project(HMOranges)
set(CMAKE_CXX_STANDARD 11)

MESSAGE(Compiler: \ ${CMAKE_CXX_COMPILER_ID} \ ${CMAKE_CXX_COMPILER_VERSION})
MESSAGE(CMake\ Version: \ ${CMAKE_VERSION})

# Fix Windows path environment variable.
if (WIN32)
    if(DEFINED ENV{OCL_ROOT})
        STRING(REGEX REPLACE "\\\\" "/" OCL_ROOT $ENV{OCL_ROOT})
    endif()
    if(DEFINED ENV{CUDA_PATH})
        STRING(REGEX REPLACE "\\\\" "/" CUDA_PATH $ENV{CUDA_PATH})
    endif()
endif(WIN32)
if(DEFINED OCL_ROOT AND DEFINED CUDA_PATH)
    MESSAGE("Both AMD and NVIDIA OpenCL SDKs have been found. Defaulting to AMD.")
endif()
if (DEFINED OCL_ROOT)
    MESSAGE("Building with AMD OpenCL SDK.")
    # Set required OpenCL CMake build options.
    set(OpenCL_INCLUDE_DIR ${OCL_ROOT}/include/)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        # 64 bits
        set(OpenCL_LIBRARY ${OCL_ROOT}/lib/x86_64/OpenCL.lib)
    elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
        # 32 bits
        set(OpenCL_LIBRARY ${OCL_ROOT}/lib/x86/OpenCL.lib)
    endif()

    # Add OpenCL stuff.
    find_package(OpenCL REQUIRED)
    #add_library(OpenCL ${OpenCL_LIBRARY})
    include_directories(${OpenCL_INCLUDE_DIR})
    #SET_TARGET_PROPERTIES(OpenCL PROPERTIES LINKER_LANGUAGE C)
elseif(DEFINED CUDA_PATH)
    MESSAGE("Building with NVIDIA OpenCL SDK (CUDA Toolkit).")
    set(OpenCL_INCLUDE_DIR ${CUDA_PATH}/include/)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        # 64 bits
        set(OpenCL_LIBRARY ${CUDA_PATH}/lib/x64/OpenCL.lib)
    elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
        # 32 bits
        set(OpenCL_LIBRARY ${CUDA_PATH}/lib/Win32/OpenCL.lib)
    endif()


    # Add OpenCL stuff.
    find_package(OpenCL REQUIRED)
    #add_library(OpenCL ${OpenCL_LIBRARY})
    include_directories(${OpenCL_INCLUDE_DIR})
    #SET_TARGET_PROPERTIES(OpenCL PROPERTIES LINKER_LANGUAGE C)
else()
    MESSAGE(OCL_ROOT\ not\ found.)
    find_package(OpenCL REQUIRED)
    set(OpenCL_LIBRARY OpenCL::OpenCL)
endif()

# Allow PROJECT_DIR to be used in C source files (useful for runtime kernel compilation)
add_compile_definitions(PROJECT_DIR="${PROJECT_SOURCE_DIR}/")

# General source files

file(GLOB_RECURSE sim_SOURCE_FILES
        "util/*.h" "util/*.c"
        "util/particleUtils/particleUtils.c" "util/particleUtils/particleUtils.h"
        "structures/*.h" "util/kernelUtils.cl"
        "tests/*/*"
        "sims/simRunner/*")

# Targets

file(GLOB_RECURSE run_tests_SOURCE_FILES
        "util/*.h" "util/*.c"
        "structures/*.h"
        tests/*)
add_executable(run_tests ${run_tests_SOURCE_FILES} tests/standalone_tests.c)
target_link_libraries(run_tests ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(run_tests m)
endif()

file(GLOB_RECURSE tgv_periodic_SOURCE_FILES sims/tgv_periodic.c)
add_executable(tgv_periodic ${tgv_periodic_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(tgv_periodic ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(tgv_periodic m)
endif()

file(GLOB_RECURSE multi_tgv_periodic_SOURCE_FILES sims/multi_tgv_periodic.c)
add_executable(multi_tgv_periodic ${multi_tgv_periodic_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(multi_tgv_periodic ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(multi_tgv_periodic m)
endif()

file(GLOB_RECURSE drag_SOURCE_FILES verification/drag/drag_sim.c)
add_executable(drag ${drag_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(drag ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(drag m)
endif()

file(GLOB_RECURSE gravity_drag_SOURCE_FILES verification/gravity_drag/gravity_drag_sim.c)
add_executable(gravity_drag ${gravity_drag_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(gravity_drag ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(gravity_drag m)
endif()

file(GLOB_RECURSE uncoupled_heat_transfer_SOURCE_FILES verification/uncoupled_heat_transfer/uc_heat_transfer_sim.c)
add_executable(uncoupled_heat_transfer ${uncoupled_heat_transfer_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(uncoupled_heat_transfer ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(uncoupled_heat_transfer m)
endif()

file(GLOB_RECURSE uncoupled_mass_transfer_SOURCE_FILES verification/uncoupled_mass_transfer/uc_mass_transfer_sim.c)
add_executable(uncoupled_mass_transfer ${uncoupled_mass_transfer_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(uncoupled_mass_transfer ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(uncoupled_mass_transfer m)
endif()

file(GLOB_RECURSE coupled_heat_mass_transfer_SOURCE_FILES verification/coupled_heat_mass_transfer/c_heat_mass_transfer_sim.c)
add_executable(coupled_heat_mass_transfer ${coupled_heat_mass_transfer_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(coupled_heat_mass_transfer ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(coupled_heat_mass_transfer m)
endif()

file(GLOB_RECURSE stokes_numbers_SOURCE_FILES analysis/stokes_numbers/stokes_numbers_sim.c)
add_executable(stokes_numbers ${stokes_numbers_SOURCE_FILES} ${sim_SOURCE_FILES})
target_link_libraries(stokes_numbers ${OpenCL_LIBRARY})
if (NOT(WIN32))
    target_link_libraries(stokes_numbers m)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -funroll-loops ") # Optimization stuff, TODO: Investigate.
#SET(CMAKE_VERBOSE_MAKEFILE ON)
