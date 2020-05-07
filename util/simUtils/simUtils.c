//
// Created by Elijah on 15/02/2018.
//

#include "simUtils.h"

int writeParticles(particle *particles, float time, char prefix[], char dir[], cl_ulong NUMPART, boolean log_vel) {
    char filename[500];

    if (!checkDirExists(dir)) {
        fprintf(stderr, "Error: Directory (%s) does not exist or cannot be accessed.\n", dir);
        return 0;
    }

    sprintf(filename, "%s%llu_%s_%i.txt", dir, NUMPART, prefix, (long) roundf(time * 1e6));
    FILE *fd = fopen(filename, "w");
    for (int i = 0; i < NUMPART; i++) {
        if (log_vel) {
            fprintf(fd, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\n", particles[i].pos.x, particles[i].pos.y, particles[i].pos.z,
                    particles[i].vel.x, particles[i].vel.y, particles[i].vel.z, particles[i].fluid_vel.x,
                    particles[i].fluid_vel.y, particles[i].fluid_vel.z, particles[i].T_d, particles[i].diameter*1000*particles[i].diameter*1000);
        } else {
            fprintf(fd, "%f,%f,%f\n", particles[i].pos.x, particles[i].pos.y, particles[i].pos.z);
        }
    }
    fclose(fd);
    return 1;
}

int writeSetupData(char prefix[], char dir[], cl_ulong NUMPART, cl_float timestep, cl_float sim_length,
                   cl_float domain_length,
                   cl_float particle_diameter, cl_float effect_diameter, cl_float particle_density,
                   cl_float fluid_viscosity) {
    char filename[500];

    if (!checkDirExists(dir)) {
        fprintf(stderr, "Error: Directory (%s) does not exist or cannot be accessed.\n", dir);
        return 0;
    }

    sprintf(filename, "%s%llu_%s_setup.txt", dir, NUMPART, prefix);
    FILE *fd = fopen(filename, "w");
    fprintf(fd, "Number of particles: %llu\n", NUMPART);
    fprintf(fd, "Timestep: %f\n", timestep);
    fprintf(fd, "Simulation length: %f\n", sim_length);
    fprintf(fd, "Domain length: %f\n", domain_length);
    fprintf(fd, "Particle diameter: %f\n", particle_diameter);
    fprintf(fd, "Effect diameter: %f\n", effect_diameter);
    fprintf(fd, "Particle density: %f\n", particle_density);
    fprintf(fd, "Fluid viscosity: %f\n", fluid_viscosity);
    fclose(fd);
    return 1;
}

int writeTime(char prefix[], char dir[], cl_ulong NUMPART, char label[]) {
    char filename[500];

    if (!checkDirExists(dir)) {
        fprintf(stderr, "Error: Directory (%s) does not exist or cannot be accessed.\n", dir);
        return 0;
    }

    sprintf(filename, "%s%llu_%s_setup.txt", dir, NUMPART, prefix);

    FILE *fd = fopen(filename, "a");

    time_t now;
    now = time(NULL);
    fprintf(fd, "%s: %s\n", label, ctime(&now));
    fclose(fd);
    return 1;
}

boolean checkDirExists(char dir[]) {
    if (strcmp(dir, "") == 0) {
        return TRUE;
    }
#if defined(_MSC_VER)
    struct stat sb;
    if (GetFileAttributes(dir) == INVALID_FILE_ATTRIBUTES) {
        return FALSE;
    } else {
        return TRUE;
    }
#elif defined(__GNUC__) || defined(__GNUG__) || defined(__MINGW_GCC_VERSION)
    // https://stackoverflow.com/a/12510903/5270376
    DIR *dir_obj = opendir(dir);
    if (dir_obj) {
        /* Directory exists. */
        closedir(dir_obj);
        return TRUE;
    } else if (ENOENT == errno) {
        /* Directory does not exist. */
        return FALSE;
    } else {
        /* opendir() failed for some other reason. */
        return FALSE;
    }
#endif
}

float getFluidViscFromStks(float diameter, float density, float flow_speed, float length, float stks) {
    return density * diameter * diameter * flow_speed / (18 * stks * length);
}

float getRestitutionFromSy(float eff_len, float k_c, float u, float mass, float sy) {
    return eff_len * sqrtf(k_c) / (sy * sqrtf(eff_len * eff_len * k_c + u * u * mass));
}

float getCohesionFromSy(float eff_len, float restitution, float u, float mass, float sy) {
    return powf(sy * restitution * u, 2) * mass / (powf(eff_len, 2) * (1 - powf(sy * restitution, 2)));
}