__kernel void test_particle_struct_alignment(__global particle *particles, __global bool *correct) {
        int gid = get_global_id(0);

        if (!(particles[gid].pos.x == 20 &&
              particles[gid].pos.y == 21 &&
              particles[gid].pos.z == 22 &&
              particles[gid].vel.x == 23 &&
              particles[gid].vel.y == 24 &&
              particles[gid].vel.z == 25 &&
              particles[gid].forces.x == 26 &&
              particles[gid].forces.y == 27 &&
              particles[gid].forces.z == 28 &&
              particles[gid].fluid_vel.x == 29 &&
              particles[gid].fluid_vel.y == 30 &&
              particles[gid].fluid_vel.z == 31 &&
              particles[gid].id == 32 &&
              particles[gid].cv_array_idx == 33 &&
              particles[gid].diameter == 34 &&
              particles[gid].effect_diameter == 35 &&
              particles[gid].density == 36 &&
              particles[gid].fluid_viscosity == 37 &&
              particles[gid].T_d == 38 &&
              particles[gid].W_V == 39 &&
              particles[gid].T_B == 40 &&
              particles[gid].L_V == 41 &&
              particles[gid].C_L == 42 &&
              particles[gid].P_atm == 43 &&
              particles[gid].R_bar == 44 &&
              particles[gid].R == 45 &&
              particles[gid].W_G == 46 &&
              particles[gid].theta_1 == 47 &&
              particles[gid].theta_2 == 48 &&
              particles[gid].Y_G == 49 &&
              particles[gid].Pr_G == 50 &&
              particles[gid].Sc_G == 51 &&
              particles[gid].f_2 == 52 &&
              particles[gid].P_G == 53 &&
              particles[gid].T_G == 54 &&
              particles[gid].H_deltaT == 55 &&
              particles[gid].m_d == 56
              )) {
            *correct = false;
        }

        particles[gid].pos.x = 60;
        particles[gid].pos.y = 61;
        particles[gid].pos.z = 62;
        particles[gid].vel.x = 63;
        particles[gid].vel.y = 64;
        particles[gid].vel.z = 65;
        particles[gid].forces.x = 66;
        particles[gid].forces.y = 67;
        particles[gid].forces.z = 68;
        particles[gid].fluid_vel.x = 69;
        particles[gid].fluid_vel.y = 70;
        particles[gid].fluid_vel.z = 71;
        particles[gid].id = 72;
        particles[gid].cv_array_idx = 73;
        particles[gid].diameter = 74;
        particles[gid].effect_diameter = 75;
        particles[gid].density = 76;
        particles[gid].fluid_viscosity = 77;
        particles[gid].T_d = 78;
        particles[gid].W_V = 79;
        particles[gid].T_B = 80;
        particles[gid].L_V = 81;
        particles[gid].C_L = 82;
        particles[gid].P_atm = 83;
        particles[gid].R_bar = 84;
        particles[gid].R = 85;
        particles[gid].W_G = 86;
        particles[gid].theta_1 = 87;
        particles[gid].theta_2 = 88;
        particles[gid].Y_G = 89;
        particles[gid].Pr_G = 90;
        particles[gid].Sc_G = 91;
        particles[gid].f_2 = 92;
        particles[gid].P_G = 93;
        particles[gid].T_G = 94;
        particles[gid].H_deltaT = 95;
        particles[gid].m_d = 96;
}
