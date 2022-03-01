#include "simulation.hpp"

using namespace cgp;

const vec3 g = { 0,0,-9.81f };
const float alpha = 0.9;// attenuation for collision
const float beta = 0.9;// attenuation for collision
const float epsilon = 0.001f;

// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters)
{
    // direct access to the variables
    grid_2D<vec3>& force = cloth.force;
    grid_2D<vec3> const& position = cloth.position;
    grid_2D<vec3> const& velocity = cloth.velocity;

    size_t const N = cloth.position.size();       // total number of vertices
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid

    // Retrive simulation parameter
    float const K = parameters.K;              // spring stifness
    float const m = parameters.mass_total / N; // mass of a particle
    float const mu = parameters.mu;            // damping coefficient
    float const	L0 = 1.0f / (N_edge - 1.0f);   // rest length between two direct neighboring particle

    // Gravity
    for (int ku = 0; ku < N_edge; ++ku)
        for (int kv = 0; kv < N_edge; ++kv)
            force(ku, kv) = m * g;

    // Drag
    for (int ku = 0; ku < N_edge; ++ku)
        for (int kv = 0; kv < N_edge; ++kv)
            force(ku, kv) += -mu * m * velocity(ku, kv);


    // TO DO: Add spring forces ...
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            const vec3& pi = position(ku, kv);
            vec3 sum(0, 0, 0);
            vec3 pij;
            int index = position.index_to_offset(ku, kv);
            // neighbor and bending springs
            for (int3 neighborInfo : cloth.neighbor_cross_map[index]) {
                float Lij = L0 * neighborInfo[2];
                pij = position(neighborInfo[0], neighborInfo[1]) - pi;
                sum += (1 - Lij / norm(pij)) * pij;
            }
            // shearing springs
            for (int3 neighborInfo : cloth.neighbor_diagonal_map[index]) {
                float Lij = sqrt(2) * L0 * neighborInfo[2];
                pij = position(neighborInfo[0], neighborInfo[1]) - pi;
                sum += (1 - Lij / norm(pij)) * pij;
            }

            // add wind force
            // F_wind = dot(normal, v_wind)*Area*m* normal* alpha
            vec3 surface_normal = cloth.normal(ku, kv);
            vec3 wind_dir = parameters.wind.direction;
            if (norm(surface_normal) >= 1e-5f && norm(wind_dir) >= 1e-5f) {
                surface_normal = normalize(surface_normal);
                wind_dir = normalize(wind_dir);
                float n_dot_wind_dir = dot(surface_normal, wind_dir);
                vec3 windforce_dir = n_dot_wind_dir > 0 ? surface_normal : -surface_normal;
                vec3 force_wind = std::abs(n_dot_wind_dir) * windforce_dir * m * parameters.wind.magnitude;
                force(ku, kv) += force_wind;
            }

            force(ku, kv) += K * sum;
        }
    }
}

void simulation_numerical_integration(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters const& parameters, float dt)
{
    int const N_edge = cloth.N_samples_edge();
    float const m = parameters.mass_total / static_cast<float>(N_edge);

    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            vec3& v = cloth.velocity(ku, kv);
            vec3& p = cloth.position(ku, kv);
            vec3 const& f = cloth.force(ku, kv);

            // Standard semi-implicit numerical integration
            v = v + dt * f / m;
            p = p + dt * v;
        }
    }

    vec3 const f = 0.05f * g;

    falling_sphere.v = (1 - 0.9f * dt) * falling_sphere.v + dt * f;
    falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
}

void simulation_apply_constraints(cloth_structure& cloth, particle_structure& falling_sphere, constraint_structure const& constraint, simulation_parameters const& parameters)
{
    // Fixed positions of the cloth
    for (auto const& it : constraint.fixed_sample) {
        position_contraint c = it.second;
        cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
    }


    // To do: apply external constraints
    // For all vertex:
    //   If vertex is below floor level ...
    //   If vertex is inside collision sphere ...
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    const float radius_offset = 0.01f; // because we just check the vertices not triangles, that will lead to vertices on sphere's surface but the middle point get inside into sphere
    float const m = parameters.mass_total / cloth.position.size(); // mass of a particle
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            vec3& p = cloth.position(ku, kv);
            
            // collide with spheres
            vec3 normal = cloth.normal(ku, kv);
            vec3 cp = p - constraint.sphere.center; // center to p
            float cp_len = norm(cp);
            float diff = cp_len - (constraint.sphere.radius + radius_offset);
            if (diff <= 0) {
                //p += cp / cp_len * diff; // move along with center to p direction -> not good, will cause oscillation
                p += normal * std::abs(dot(cp / cp_len * diff, normal)); // move along with center to p direction
            }


            // YUAN test, cloth collides with falling_sphere
            const float mu = 0.8f;
            float detection = norm(p - falling_sphere.p);
            vec3 uij = (p - falling_sphere.p) / detection; // direction from pi(cloth particle) to pj(falling_sphere)
            if (detection <= (falling_sphere.r + radius_offset)) {
                vec3& vi = cloth.velocity(ku, kv);
                vec3 vij = falling_sphere.v - vi;
                if (norm(vij) <= epsilon) {
                    // static contact
                    vi *= mu;
                    falling_sphere.v *= mu;
                }
                else {
                    // collide
                    // impluse calculation
                    vec3 u = -uij; // from j to i
                    float j = 2 * (m * falling_sphere.m) / (m + falling_sphere.m) * dot(vij, u);
                    vec3 vi_perpendicular = dot(vi, u) * u;
                    vec3 vi_parallel = vi - vi_perpendicular;
                    vec3 vj_perpendicular = dot(falling_sphere.v, u) * u;
                    vec3 vj_parallel = falling_sphere.v - vj_perpendicular;
                    // apply impluse in perpendicular direction (perpendicular to seperating plane / along with u direction)
                    vec3 dv = 2 / (m + falling_sphere.m) * dot(vj_perpendicular - vi_perpendicular, u) * u;
                    vi_perpendicular += falling_sphere.m * dv;
                    vj_perpendicular -= m * dv;
                    //pi.v = alpha * vi_parallel + beta * vi_perpendicular;
                    //pj.v = alpha * vj_parallel + beta * vj_perpendicular;
                    vi = alpha * vi_parallel + beta * vi_perpendicular;
                    falling_sphere.v = alpha * vj_parallel + beta * vj_perpendicular;
                }

                // correct position with simplest method
                vec3 d = (falling_sphere.r - detection) / 2 * uij;
                p -= d;
                falling_sphere.p += d;
            }

            // collide with ground, add epsilon to avoid rendering cloth and plane alternately when their triangle overlap
            p.z = std::max(constraint.ground_z + epsilon, p.z);
        }
    }
}

bool simulation_detect_divergence(cloth_structure const& cloth)
{
    bool simulation_diverged = false;
    const size_t N = cloth.position.size();
    for (size_t k = 0; simulation_diverged == false && k < N; ++k)
    {
        const float f = norm(cloth.force.data.at_unsafe(k));
        const vec3& p = cloth.position.data.at_unsafe(k);

        if (std::isnan(f)) // detect NaN in force
        {
            std::cout << "\n **** NaN detected in forces" << std::endl;
            simulation_diverged = true;
        }

        if (f > 600.0f) // detect strong force magnitude
        {
            std::cout << "\n **** Warning : Strong force magnitude detected " << f << " at vertex " << k << " ****" << std::endl;
            simulation_diverged = true;
        }

        if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) // detect NaN in position
        {
            std::cout << "\n **** NaN detected in positions" << std::endl;
            simulation_diverged = true;
        }
    }

    return simulation_diverged;
}