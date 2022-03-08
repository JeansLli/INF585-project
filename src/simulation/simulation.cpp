#include "simulation.hpp"
#include <math.h>

using namespace cgp;

const vec3 g = { 0,0,-9.81f };
const float radius_offset = 0.025f; // because we just check the vertices not triangles, that will lead to vertices on sphere's surface but the middle point get inside into sphere
//const float alpha = 0.9;// attenuation for collision
//const float beta = 0.9;// attenuation for collision
const float epsilon = 0.01f;
float sphere_first_projection; // sphere center project on surface_normal for first contact with cloth

//Jingyi's work 2022-3-3
void collision_sphere_sphere(particle_structure& sphere, cloth_structure& cloth, int ku, int kv, float radius_offset)
{
    vec3& p1 = sphere.p;
    vec3& p2 = cloth.position(ku, kv);
    vec3& v1 = sphere.v;
    vec3& v2 = cloth.velocity(ku, kv);
    float const r1 = sphere.r;
    float const r2 = radius_offset;

    float const epsilon = 1e-5f;
    float const alpha = 0.95f;

    vec3 const p12 = p1 - p2;
    float const d12 = norm(p12);

    if (d12 < r1 + r2)
    {
        vec3 const u12 = p12 / d12;
        float const collision_depth = r1 + r2 - d12;

        p1 += (collision_depth / 2.0f + epsilon) * u12;
        p2 -= (collision_depth / 2.0f + epsilon) * u12;

        if (norm(v1 - v2) > 0.2f) {
            float const j = dot(v1 - v2, u12);
            v1 = v1 - alpha * j * u12;
            v2 = v2 + alpha * j * u12;
        }
        else // Contact
        {
            v1 = v1 * 0.83f;
            v2 = v2 * 0.83f;
        }
    }

    /* This's Putian's work for two spheres collision
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
    */
}


//Jingyi's work done 


void simulation_collision_detection(cloth_structure &cloth,particle_structure& falling_sphere, simulation_parameters& parameters){
    if (!parameters.need_collision_detect)
        return; // when sphere is in the seconde phrase of bouncing, no need to and should not to detect collision cause if do so, if will connect cloth again
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            if (cloth.contact_info[index].is_contact)
                continue;
            vec3& p = cloth.position(ku, kv);
            vec3 cf = p - falling_sphere.p;
            float cf_len = norm(cf);
            if(cf_len <= (falling_sphere.r + radius_offset)){
                if (!parameters.is_connecting) {
                    sphere_first_projection = dot(falling_sphere.p, cloth.surface_normal);
                    falling_sphere.v *= 0.7f; // regression, velocity loss in collision
                }
                parameters.is_connecting = true;
                parameters.is_extending = true;
                cloth.contact_info[index].is_contact = true;
                cloth.contact_info[index].offset = cgp::normalize(cf) * (falling_sphere.r + radius_offset);
            }
        }
    }         
}

void fall_sphere_update(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters& parameters, float dt){
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    vec3 const g_scale = 0.05f * g;
    if(!parameters.is_connecting){
        //std::cout << "fall freely"<<std::endl;
        falling_sphere.v = falling_sphere.v + dt * g_scale;
        falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
        //std::cout << "no extension, v=" << falling_sphere.v << std::endl;
        float v_dot_sn = dot(falling_sphere.v, cloth.surface_normal);
        parameters.need_collision_detect = v_dot_sn <= 0 && std::abs(v_dot_sn) >= epsilon; // if sphere is moving forward to the cloth plane
    }
    else {
        // TODO: extend fall_sphere to multiple sphere, handle them together, put all corresponding data into one single data structure
        const float K = 50.0f;
        const float sphere_projection = dot(falling_sphere.p, cloth.surface_normal);
        const float offset = sphere_first_projection - sphere_projection;
        const cgp::vec3 f = K * std::max(offset, 0.0f) * cloth.surface_normal;
        const cgp::vec3 a = f / falling_sphere.m + g_scale;

        float v_projection = dot(falling_sphere.v, cloth.surface_normal);
        float a_projection = dot(a, cloth.surface_normal);
        if (std::abs(a_projection) <= epsilon && v_projection >= 0 && v_projection <= epsilon) {
            //static contact
            return;
        }

        if (parameters.is_extending) {
            if (v_projection >= 0) {
                parameters.is_extending = false;
                parameters.need_collision_detect = false; // in bouncing-back, no need to collision detection anymore

                //// mirror reflect
                //vec3 v_n = v_n_length * cloth.surface_normal;
                //vec3 v_parallel = falling_sphere.v - v_n;
                //falling_sphere.v = v_parallel - v_n;

                //detach all constraints particles
                for (int k = 0; k < cloth.contact_info.size(); k++) {
                    cloth.contact_info[k].is_contact = false;
                }
            }
            else {
                // using string force
                falling_sphere.v = falling_sphere.v + a * dt;
                falling_sphere.p = falling_sphere.p + falling_sphere.v * dt;
            }
        }
        else {
            if (offset <= 0) {
                parameters.is_connecting = false; // finish bouncing
                parameters.need_collision_detect = false; // in bouncing-back, no need to collision detection anymore
            }
            else {
                // using string force
                const float coef_loss = 0.98f;  // regression, must be regression, just for avoid oscillation
                falling_sphere.v = coef_loss * falling_sphere.v + a * dt;
                falling_sphere.p = falling_sphere.p + falling_sphere.v * dt;
            }
        }
    }
}

void update_cloth_constraints(cloth_structure& cloth, particle_structure& falling_sphere){
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            if(cloth.contact_info[index].is_contact){
                vec3& p = cloth.position(ku, kv);
                p = falling_sphere.p + cloth.contact_info[index].offset;
            }
        }
    }

}


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


    // Add spring forces ...
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
    
    //Jingyi's work 2022-3-3
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            if (cloth.contact_info[index].is_contact)
                continue;
            vec3& v = cloth.velocity(ku, kv);
            vec3& p = cloth.position(ku, kv);
            vec3 const& f = cloth.force(ku, kv);
            // Standard semi-implicit numerical integration
            v = v + dt * f / m;
            p = p + dt * v;
        }
    }
//Jingyi's work done
}


void simulation_apply_constraints(cloth_structure& cloth, particle_structure& falling_sphere, constraint_structure const& constraint, simulation_parameters const& parameters)
{
    // Fixed positions of the cloth
    for (auto const& it : constraint.fixed_sample) {
        position_contraint c = it.second;
        cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
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