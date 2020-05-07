#ifndef DEFAULT_PARTICLE_GRAVITY
float3 get_gravity(particle p, float delta_t) {
    return (float3) {0, -9.81, 0};
}
#endif

#ifndef DEFAULT_PARTICLE_FLUID_VEL
float3 get_vel_fluid(particle p, float delta_t) {
    return (float3) {0, 0, 0};
}
#endif

float get_tau(particle p) {
    return p.density * p.diameter * p.diameter / (18.0 * p.fluid_viscosity);
}

float get_u_s(particle p, float delta_t){
    float3 fluid_vel = get_vel_fluid(p, delta_t);
    float3 rel_vel = {fluid_vel.x-p.vel.x, fluid_vel.y-p.vel.y, fluid_vel.z-p.vel.z};
    return sqrt(rel_vel.x * rel_vel.x + rel_vel.y * rel_vel.y + rel_vel.z * rel_vel.z);
}

double get_Re(particle p, float delta_t) {
    double u_s = get_u_s(p, delta_t);
    double mu_G = p.fluid_viscosity;
    double rho_G = p.rho_G;
    return (rho_G*u_s*p.diameter)/(mu_G);
}

double get_Sh(particle p, float delta_t, int fixed_Re){
    double Re_d;
    double power_one_three = 1.0/3.0;
    double power_one_two = 1.0/2.0;
    if (fixed_Re == 0) {
        Re_d = get_Re(p, delta_t);
    } if (fixed_Re == 1) {
        Re_d = p.Reynolds_d;
    }
    return 2 + (0.552*pow(Re_d, power_one_two)*pow(p.Sc_G, power_one_three));
}

double get_Nu(particle p, float delta_t, int fixed_Re){
    double Re_d;
    double power_one_three = 1.0/3.0;
    double power_one_two = 1.0/2.0;
    if (fixed_Re == 0) {
        Re_d = get_Re(p, delta_t);
    } if (fixed_Re == 1) {
        Re_d = p.Reynolds_d;
    }
    return 2 + (0.552*pow(Re_d, power_one_two)*pow(p.Pr_G, power_one_three));
}

double get_mass_fractions(particle p){
    double x_s_eq = (p.P_atm/p.P_G) * exp(((p.L_V)/(p.R_bar/p.W_V))*((1/p.T_B)-(1/p.T_d)));
    double Y_s_eq = (x_s_eq) / (x_s_eq+(1-x_s_eq)*p.theta_2);
    double B_m_eq = (Y_s_eq-p.Y_G) / (1-Y_s_eq);
    double H_M = log(1.0 + B_m_eq);
    return H_M;
}

double get_m_d_dot(particle p, float delta_t, int fixed_Re) {
    double tau = get_tau(p);
    double Sh = get_Sh(p, delta_t, fixed_Re);
    double H_M = get_mass_fractions(p);
    return -((((Sh)/(3*p.Sc_G))*(p.m_d/tau))*H_M);
}

float3 get_accel(particle p, float delta_t) {
    float tau = get_tau(p);
    float3 non_drag_a = get_gravity(p, delta_t);
    return (get_vel_fluid(p, delta_t) - p.vel + tau * non_drag_a) / (tau + delta_t);
}

double iterate_temp(particle p, float delta_t, int coupled, int analytic, int fixed_Re) {
    double temp_return;
    double tau = get_tau(p);
    double Nu = get_Nu(p, delta_t, fixed_Re);
    double m_d_dot = get_m_d_dot(p, delta_t, fixed_Re);
    if (coupled == 1 && analytic == 1) {
        temp_return = p.T_G - ((p.T_G - p.T_d)*exp( -((((p.f_2*Nu)/(3.0*p.Pr_G))*(p.theta_1/tau))*delta_t)));
    } if (coupled == 1 && analytic == 2) {
        temp_return = ((p.T_d + delta_t*((((p.f_2*Nu)/(3.0*p.Pr_G))*(p.theta_1/tau)*(p.T_G))))/(1.0+delta_t*(((p.f_2*Nu)/(3.0*p.Pr_G))* (p.theta_1/tau))));
    } if (coupled == 2 && analytic == 2){
        temp_return = (p.T_d + delta_t*((((p.f_2*Nu)/(3.0*p.Pr_G)) * (p.theta_1/tau)*(p.T_G)) + ((p.L_V*m_d_dot)/(p.C_L*p.m_d))-p.H_deltaT)) / (1.0+delta_t*(((p.f_2*Nu)/(3.0*p.Pr_G)) * (p.theta_1/tau)));
    }
    return (temp_return);
}

// uc ((p.T_d + delta_t*((((p.f_2*Nu)/(3*p.Pr_G))*(p.theta_1/tau)*(p.T_G))))/ (1+delta_t*(((p.f_2*Nu)/(3*p.Pr_G))* (p.theta_1/tau))))

// c (p.T_d + delta_t*((((p.f_2*Nu)/(3*p.Pr_G)) * (p.theta_1/tau)*(p.T_G)) + ((p.L_V*m_d_dot)/(p.C_L*m_d))-p.H_deltaT)) / (1+delta_t*(((p.f_2*Nu)/(3*p.Pr_G)) * (p.theta_1/tau)))

double iterate_mass(particle p, float delta_t, int coupled, int analytic, int fixed_Re){
    double mass_return;
    double mass_for_pow;
    double density_for_pow;
    double power_3_2 = 3.0/2.0;
    double power_2_3 = 2.0/3.0;
    double power_1_3 = 1.0/3.0;
    double tau = get_tau(p);
    double Sh = get_Sh(p, delta_t, fixed_Re);
    double H_M = get_mass_fractions(p);
    if (analytic == 2){
        mass_return = (p.m_d)/(1.0 + (delta_t*(((Sh)/(3.0*p.Sc_G))*(H_M/tau))));
    } else{
        density_for_pow = 6.0 / p.density;
        mass_for_pow = pow(p.m_d, 2.0/3.0); - (delta_t*(((2.0*Sh*p.fluid_viscosity*H_M)/(p.Sc_G*3.0)*pow(density_for_pow, 1.0/3.0)*pow(M_PI, 2.0/3.0))));
        mass_return = pow(mass_for_pow, 3.0/2.0);
    }
    return (mass_return);

}

float3 iterate_velocity(particle p, float delta_t) {
    float3 next_vel = p.vel + delta_t * get_accel(p, delta_t);
    return next_vel;
}

float3 iterate_position(particle p, float delta_t, float3 next_vel) {
    return p.pos + delta_t * (next_vel + p.vel) / 2;
}

/* Kernel to iterate particles. */
__kernel void iterate_particle(__global particle *particles, float delta_t, int int_periodic, float domain_length,
        int coupled, int analytic, int fixed_Re) {
    bool periodic = (bool) int_periodic;
    int gid = get_global_id(0);
    if (particles[gid].density == -1) {
        return; // -1 is used to denote infinite density.
    }
    double D_0 = (pow(1.1, 0.5)/1000);
    double md_0 = 997 * M_PI * (D_0 * D_0 * D_0) / 6;

    if (coupled == 0 && analytic == 0){
        if (particles[gid].m_d/md_0 <= 0.0001){
            particles[gid].pos = 0;
            particles[gid].vel = 0;
            particles[gid].diameter = 0;
            particles[gid].T_d = 0;
        } else{
            float next_mass = iterate_mass(particles[gid], delta_t, coupled,  analytic, fixed_Re);
            particles[gid].m_d = next_mass;
            double mass_for_pow = (next_mass * 6.0) / (particles[gid].density * M_PI);
            particles[gid].diameter = pow(mass_for_pow, 1.0/3.0);
        }
    }

    if (coupled == 1 && analytic == 1){
        double next_temp = iterate_temp(particles[gid], delta_t, coupled, analytic, fixed_Re);
        particles[gid].T_d = next_temp;
    }

    if (coupled == 0 && analytic == 2){
        if (particles[gid].diameter == 0){
            particles[gid].pos = 0;
            particles[gid].vel = 0;
            particles[gid].diameter = 0;
            particles[gid].T_d = 0;
        } else{
            double next_mass = iterate_mass(particles[gid], delta_t, coupled, analytic, fixed_Re);
            particles[gid].m_d = next_mass;
            particles[gid].diameter = pow((next_mass*6)/(particles[gid].density*M_PI),1.0/3.0);
        }
    }

    if (coupled == 1 && analytic == 2){
        if (particles[gid].T_d == particles[gid].T_G){
            particles[gid].pos = 0;
            particles[gid].vel = 0;
            particles[gid].diameter = 0;
            particles[gid].T_d = 0;
        } else{
            double next_temp = iterate_temp(particles[gid], delta_t, coupled, analytic, fixed_Re);
            particles[gid].T_d = next_temp;
        }
    }

    if (coupled == 2 && analytic == 2 && fixed_Re == 1){
        if (particles[gid].m_d/particles[gid].initial_mass <= 0.001){
            particles[gid].pos = 0;
            particles[gid].vel = 0;
            particles[gid].diameter = 0;
            particles[gid].T_d = 0;
        } else { //if (particles[gid].diameter/0.001048 > 0.001){ //(particles[gid].mass/particles[gid].initial_mass > 0.001
            double next_mass = iterate_mass(particles[gid], delta_t, coupled, analytic, fixed_Re);
            double next_temp = iterate_temp(particles[gid], delta_t, coupled, analytic, fixed_Re);
            particles[gid].m_d = next_mass;
            double mass_for_pow = (next_mass * 6.0) / (particles[gid].density * M_PI);
            particles[gid].diameter = pow(mass_for_pow, 1.0/3.0);
            particles[gid].T_d = next_temp;
        }
    }

    if (coupled == 2 && analytic == 2 && fixed_Re == 0){
        if (particles[gid].m_d/particles[gid].initial_mass <= 0.001){
            particles[gid].pos = 0;
            particles[gid].vel = 0;
            particles[gid].diameter = 0;
            particles[gid].T_d = 0;
        } else{ //if (particles[gid].diameter/0.001048 > 0.001){ //(particles[gid].mass/particles[gid].initial_mass > 0.001
            float3 next_vel = iterate_velocity(particles[gid], delta_t);
            float3 next_pos = iterate_position(particles[gid], delta_t, next_vel);
            double next_mass = iterate_mass(particles[gid], delta_t, coupled, analytic, fixed_Re);
            double next_temp = iterate_temp(particles[gid], delta_t, coupled, analytic, fixed_Re);

            particles[gid].pos = next_pos;
            particles[gid].vel = next_vel;
            particles[gid].m_d = next_mass;
            double b = 1.0/3.0;
            particles[gid].diameter = pow((next_mass*6.0)/(particles[gid].density*M_PI), b);
            particles[gid].T_d = next_temp;

            particles[gid].fluid_vel = get_vel_fluid(particles[gid], delta_t);
            particles[gid].forces = (float3) {0, 0, 0};
        }
    }

    if (periodic) {
        if (particles[gid].pos.x > domain_length / 2) {
            particles[gid].pos.x -= domain_length;
        } else if (particles[gid].pos.x < - domain_length / 2) {
            particles[gid].pos.x += domain_length;
        }

        if (particles[gid].pos.y > domain_length / 2) {
            particles[gid].pos.y -= domain_length;
        } else if (particles[gid].pos.y < - domain_length / 2) {
            particles[gid].pos.y += domain_length;
        }

        if (particles[gid].pos.z > domain_length / 2) {
            particles[gid].pos.z -= domain_length;
        } else if (particles[gid].pos.z < - domain_length / 2) {
            particles[gid].pos.z += domain_length;
        }
    }
}