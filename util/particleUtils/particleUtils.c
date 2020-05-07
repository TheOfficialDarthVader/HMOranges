//
// Created by Elijah on 18/12/2017.
//

#include "particleUtils.h"

void printParticle(particle *part) {
    particle p = *part;
    printf("%i\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n", (int) p.id, p.diameter, p.density,
           p.fluid_viscosity, p.pos.x, p.pos.y, p.pos.z, p.vel.x, p.vel.y, p.vel.z, p.forces.x, p.forces.y, p.forces.z);
}

float get_particle_mass_from_values(float density, float diameter) {
    return (float) density * PI * powf(diameter, 3) / 6;
}

float get_particle_mass(particle *p) {
    return (float) ((*p).density * PI * powf((*p).diameter, 3) / 6);
}

float get_tau(particle *p) {
    return (*p).density * (*p).diameter * (*p).diameter / (18 * (*p).fluid_viscosity);
}

boolean checkPositions(particle *particles, cl_ulong NUMPARTS, cl_float domain_length) {
    boolean correct = TRUE;
    for (cl_ulong i = 0; i < NUMPARTS; i++) {
        if (fabs(particles[i].pos.x) > domain_length / 2
            || fabs(particles[i].pos.y) > domain_length / 2
            || fabs(particles[i].pos.z) > domain_length / 2) {
            correct = FALSE;
        }
    }
    return correct;
}

void initializeMonodisperseParticles(particle *particles, cl_ulong NUMPART, float density, float fluid_viscosity,
                                     float diameter, float effect_diameter, float C_p_G, float C_L, float P_atm, float rho_G, float R_bar,
                                     float R, float W_G, float W_V, float Y_G, float T_d, float T_B, float T_G, cl_float3 positions[],
                                     cl_float3 velocities[]) {
    for (cl_ulong i = 0; i < NUMPART; i++) {
        particles[i].id = i;
        particles[i].density = density;
        particles[i].diameter = diameter;
        particles[i].effect_diameter = effect_diameter;
        particles[i].pos = positions[i];
        if (velocities != NULL) {
            particles[i].vel = velocities[i];
        } else {
            particles[i].vel = (cl_float3) {0.0, 0.0, 0.0};
        }
        particles[i].forces = (cl_float3) {0.0, 0.0, 0.0};
        particles[i].fluid_vel = (cl_float3) {0.0, 0.0, 0.0};
        particles[i].m_d = particles[i].density * M_PI * (particles[i].diameter * particles[i].diameter * particles[i].diameter) / 6;
        particles[i].initial_mass = particles[i].density * M_PI * (particles[i].diameter * particles[i].diameter * particles[i].diameter) / 6;

        //Create reference temperature//
        double T_R = (137 * pow((T_B/373.15), 0.68) * log10(T_G)) - 45;

        //droplet/liquid properties
        particles[i].C_L = C_L;
        particles[i].T_d = T_d;
        particles[i].T_B = T_B;
        particles[i].Reynolds_d = 0;

        //vapour properties
        particles[i].L_V = (2.257*(pow(10,6)) + (2.595*(pow(10,3)) * (373.15-T_G)));
        particles[i].W_V = W_V;

        //gas properties
        particles[i].P_atm = P_atm;
        particles[i].R_bar = R_bar;
        particles[i].R = R;
        particles[i].Y_G = Y_G;
        particles[i].T_G = T_G;
        particles[i].rho_G = rho_G;
        particles[i].P_G = rho_G*R*T_R;
        particles[i].fluid_viscosity = fluid_viscosity;// (6.109*pow(10,-6)) + (4.604*pow(10,-8)*T_R) - (1.05*pow(10, -11)*T_R*T_R);

        //other properties
        particles[i].H_deltaT = 0;
        particles[i].f_2 = 1;
        particles[i].theta_1 = C_p_G/C_L;
        particles[i].theta_2 = W_G/W_V;

        //non-dimensional numbers
        particles[i].Pr_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R);
        particles[i].Sc_G = 0.815 - (4.958*pow(10,-4)*T_R) + (4.514*pow(10,-7)*T_R*T_R); // =Pr_G
    }
}

float createCubePositions(cl_float3 *positions, cl_ulong NUMPART, float particle_diameter, float spacing_factor, cl_float3 center) {
    cl_ulong cubert_NUMPART = (cl_ulong) ceil(pow(NUMPART, 0.334));
    cl_ulong pos_len = 0;
    srand(0);

    float cube_length = spacing_factor * cubert_NUMPART * particle_diameter;
    for (int x = 0; x < cubert_NUMPART; x++) {
        for (int y = 0; y < cubert_NUMPART; y++) {
            for (int z = 0; z < cubert_NUMPART; z++) {
                if (pos_len < NUMPART) {
                    float rand_offset = 0.005 * (float) rand() / (float) (RAND_MAX);
                    cl_float xf = cube_length
                                  * (-0.5 + rand_offset + ((float) x / cubert_NUMPART)) + center.x;
                    cl_float yf = cube_length
                                  * (-0.5 + rand_offset + ((float) y / cubert_NUMPART)) + center.y;
                    cl_float zf = cube_length
                                  * (-0.5 + rand_offset + ((float) z / cubert_NUMPART)) + center.z;
                    positions[pos_len] = (cl_float3) {xf, yf, zf};
                }
                pos_len++;
            }
        }
    }

    return (float) cube_length;
}

void createNormalDistVelocities(cl_float3 *velocities, cl_ulong NUMPART, float mean, float std_dev) {
    srand(0);

    for (int i = 0; i < NUMPART; i++) {
        // Select random speed.
        // Box-Muller transformation from uniform to normal.
        float u = (float) (rand() + 1) / (float) (RAND_MAX + 1);
        float v = (float) (rand() + 1) / (float) (RAND_MAX + 1);
        float speed = (float) (std_dev * sqrtf(-2 * logf(u)) * cos(2 * PI * v)) + mean;

        // Select random direction.
        float pitch = (float) rand() / (float) (RAND_MAX) * PI; // From -ve z (0) to +ve z (pi).
        float yaw = (float) rand() / (float) (RAND_MAX) * 2 * PI; // +ve x = 0, -ve y = pi/2, -ve x = pi, +ve y = 3pi/4.

        float x = speed * cosf(pitch) * cosf(yaw);
        float y = speed * cosf(pitch) * sinf(yaw);
        float z = speed * sinf(pitch);

        velocities[i] = (cl_float3) {x, y, z};
    }
}