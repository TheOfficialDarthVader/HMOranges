typedef struct {
    float3 pos;
    float3 vel;
    float3 forces;
    float3 fluid_vel;
    ulong id;
    ulong cv_array_idx;
    double diameter;
    float effect_diameter;
    double density;
    double fluid_viscosity;
    double T_d;
    double W_V;
    double T_B;
    double L_V;
    double C_L;
    double P_atm;
    double R_bar;
    double R;
    double W_G;
    double theta_1;
    double theta_2;
    double Y_G;
    double Pr_G;
    double Sc_G;
    double f_2;
    double P_G;
    double T_G;
    double rho_G;
    double H_deltaT;
    double m_d;
    double Reynolds_d;
    double initial_mass;
} __attribute__ ((aligned (512))) particle;

float get_particle_mass(particle p) {
    if (p.density == -1) {
        return -1;
    } else {
        return (float) p.density * M_PI_F * pow(p.diameter, 3) / 6;
    }
}

