#define DEFAULT_PARTICLE_GRAVITY 0

float3 get_gravity(particle p, float delta_t) {
    return (float3) {0, -9.81, 0};
}
