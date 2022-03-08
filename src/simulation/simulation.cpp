#include "simulation.hpp"
#include <math.h>

using namespace cgp;

const vec3 g = { 0,0,-9.81f };
const float radius_offset = 0.015f; // because we just check the vertices not triangles, that will lead to vertices on sphere's surface but the middle point get inside into sphere
const float epsilon = 0.01f;


void collision_sphere_sphere(std::vector<sphere_structure>& spheres)
{
    const float mu = 0.8f;
    const float alpha = 0.9f;// attenuation for collision
    const float beta = 0.9f;// attenuation for collision
    const float e = 0.001f; // relative speed detection
    // particles collide with themselves
    for (size_t i = 0; i < spheres.size(); i++)
    {
        sphere_structure& pi = spheres[i];
        if (!pi.can_move)
            continue;
        for (size_t j = 0; j < spheres.size(); j++)
        {
            if (i == j) {
                continue; // no need to handle itself
            }
            sphere_structure& pj = spheres[j];
            if (!pj.can_move)
                continue;
            float detection = norm(pj.p - pi.p);
            vec3 uij = (pj.p - pi.p) / detection; // direction from pi to pj
            if (detection <= (pi.r + pj.r)) {
                vec3 vij = pj.v - pi.v;
                if (norm(vij) <= e) {
                    // static contact
                    pi.v *= mu;
                    pj.v *= mu;
                }
                else {
                    // collide
                    // impluse calculation
                    vec3 u = -uij; // from j to i
                    vec3 vi_perpendicular = dot(pi.v, u) * u;
                    vec3 vi_parallel = pi.v - vi_perpendicular;
                    vec3 vj_perpendicular = dot(pj.v, u) * u;
                    vec3 vj_parallel = pj.v - vj_perpendicular;
                    // apply impluse in perpendicular direction (perpendicular to seperating plane / along with u direction)
                    vec3 dv = 2 / (pi.m + pj.m) * dot(vj_perpendicular - vi_perpendicular, u) * u;
                    vi_perpendicular += pj.m * dv;
                    vj_perpendicular -= pi.m * dv;
                    pi.v = alpha * vi_parallel + beta * vi_perpendicular;
                    pj.v = alpha * vj_parallel + beta * vj_perpendicular;
                }

                // correct position with simplest method
                vec3 d = (pi.r + pj.r - detection) / 2 * uij;
                pi.p -= d;
                pj.p += d;
            }
        }
    }
}

void collision_sphere_cloth(cloth_structure& cloth, std::vector<sphere_structure>& spheres) {
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            cloth_contact_info& info = cloth.contact_info[index];
            if (info.connecting_spheres.size() > 0)
                // actually one particle of cloth is only contacted by one sphere at the same time
                continue;
            for (size_t i = 0; i < spheres.size(); i++) {
                sphere_structure& sphere = spheres[i];
                if (!sphere.can_move || !sphere.detect_cloth_collision)
                    // when sphere is in the seconde phrase of bouncing, no need to and should not to detect collision cause if do so, if will connect cloth again
                    continue;
                if (std::find(info.connecting_spheres.begin(), info.connecting_spheres.end(), i) != info.connecting_spheres.end())
                    // it has been handled
                    continue;
                vec3& p = cloth.position(ku, kv);
                vec3 cf = p - sphere.p;
                float cf_len = norm(cf);
                if (cf_len <= (sphere.r + radius_offset)) {
                    const float coef = 0.7f; // regression, velocity loss in collision
                    sphere.v *= coef; // regression, velocity loss in collision
                    info.offset = cgp::normalize(cf) * (sphere.r + radius_offset);
                    info.connecting_spheres.push_back(i);
                    if (sphere.connecting_particles.size() <= 0)
                        sphere.p_at_contact = sphere.p; // only record at first contact
                    sphere.connecting_particles.push_back(index);
                }
            }
        }
    }
}

void collision_spheres_box(bounding_box const& box, std::vector<sphere_structure>& spheres) {
    // spheres collide with 5 planes (not including top face)
    const float alpha = 0.9f;// attenuation for collision
    const float beta = 0.9f;// attenuation for collision
    const float half_width = box.get_width() / 2.0f;
    for (size_t k = 0; k < spheres.size(); ++k)
    {
        sphere_structure& sphere = spheres[k];
        if (!sphere.can_move)
            continue;
        for (size_t i = 0; i < box.fns.size(); i++) {
            const vec3 normal = box.get_surface_normal(i);
            const vec3 point = box.get_point_on_surface(i);
            const vec3 right = box.get_right_vector(i);
            const vec3 up = box.get_up_vector(i);

            // first, if the distance between center of sphere and face (abs(center projection on normal direction)), then sphere has possibility to collide with this face
            const vec3 offset = sphere.p - point;
            const float detection = dot(offset, normal);
            if (std::abs(detection) < sphere.r) {
                // next, check center projection on right/up direction, to see whether is inside range(width)
                const float range = half_width + sphere.r;
                const float right_projection = dot(offset, right);
                if (std::abs(right_projection) >= (range))
                    continue;
                const float up_projection = dot(offset, up);
                if (std::abs(up_projection) >= (range))
                    continue;

                // collide
                vec3 dir = detection > 0 ? normal : -normal;
                vec3 v_perpendicular = dot(sphere.v, dir) * dir;
                vec3 v_parallel = sphere.v - v_perpendicular;
                sphere.v = alpha * v_parallel - beta * v_perpendicular;

                // correct position with simplest method
                float d = sphere.r - detection;
                sphere.p += d * dir;
            }

        }
    }
}

void simulation_collision_detection(cloth_structure &cloth, std::vector<sphere_structure>& spheres, bounding_box const& box){
    // collision between spheres 
    collision_sphere_sphere(spheres);

    // collision between spheres and box
    collision_spheres_box(box, spheres);

    // collision between spheres and cloth
    collision_sphere_cloth(cloth, spheres);
}

void simulation_update_spheres(cloth_structure& cloth, std::vector<sphere_structure>& spheres, float dt){
    // update spheres position, velocity. collision b/t spheres. contact with cloth
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    vec3 const g_scale = 0.05f * g;

    const bool cloth_sphere_connecting = false;
    for (size_t i = 0; i < spheres.size(); i++) {
        sphere_structure& sphere = spheres[i];
        if (!sphere.can_move)
            continue;
        const bool is_connecting = sphere.connecting_particles.size() > 0;
        if (!is_connecting) {
            // following gravity
            sphere.v = sphere.v + dt * g_scale;
            sphere.p = sphere.p + dt * sphere.v;
            float v_projection = dot(sphere.v, vec3(0, 0, 1)); // projected on gravity direction
            sphere.detect_cloth_collision = v_projection <= 0 && std::abs(v_projection) >= epsilon; // if sphere is moving forward to the cloth plane
        }
        else {
            sphere.detect_cloth_collision = false; // in bouncing-back, no need to collision detection anymore

            vec3 surface_normal(0, 0, 0);
            for (auto index : sphere.connecting_particles) {
                int2 uv = cloth.normal.offset_to_index(index);
                surface_normal += cloth.normal(uv.x, uv.y);
            }
            surface_normal = normalize(surface_normal);
            const float K = 50.0f;
            const float p_projection_now = dot(sphere.p, surface_normal);
            const float p_projected_origin = dot(sphere.p_at_contact, surface_normal);
            const float offset = p_projected_origin - p_projection_now;
            const cgp::vec3 f = K * std::max(offset, 0.0f) * surface_normal;
            const cgp::vec3 a = f / sphere.m + g_scale;
            float v_projection = dot(sphere.v, surface_normal);
            float a_projection = dot(a, surface_normal);

            bool condition_a = a_projection <= 0 && norm(sphere.v) <= epsilon; // key condition to stop oscillation
            bool condition_b = std::abs(a_projection) <= epsilon && v_projection >= 0 && v_projection <= epsilon;
            if (condition_a || condition_b) {
                //static contact
                continue;
            }

            // TODO: optimize below codes
            const bool is_extending = v_projection < 0 && std::abs(v_projection) >= epsilon;
            if (is_extending) {
                // using string force
                sphere.v = sphere.v + a * dt;
                sphere.p = sphere.p + sphere.v * dt;
            }
            else {
                if (offset <= 0) {
                    //detach all constraints particles
                    for (size_t index : sphere.connecting_particles) {
                        cloth_contact_info& info = cloth.contact_info[index];
                        // remove sphere index
                        auto remove_it = std::remove_if(info.connecting_spheres.begin(), info.connecting_spheres.end(), [&](size_t e) { return e == i; });
                        for (auto it = remove_it; it < info.connecting_spheres.end(); it++) {
                            info.connecting_spheres.erase(it);
                        }
                    }
                    sphere.connecting_particles.clear();
                }
                else {
                    // using string force
                    const float coef_loss = 0.98f;  // regression, must be regression, just for avoid oscillation
                    sphere.v = coef_loss * sphere.v + a * dt;
                    sphere.p = sphere.p + sphere.v * dt;
                }
            }
        }
    }
}

void update_cloth_constraints(cloth_structure& cloth, std::vector<sphere_structure>& spheres){
    for (int i = 0; i < spheres.size(); i++) {
        const sphere_structure& sphere = spheres[i];
        for (int j = 0; j < sphere.connecting_particles.size(); j++) {
            size_t index = sphere.connecting_particles[j];
            int2 uv = cloth.position.offset_to_index(index);
            vec3& p = cloth.position(uv.x, uv.y);
            p = sphere.p + cloth.contact_info[index].offset;
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

void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt)
{
    int const N_edge = cloth.N_samples_edge();
    float const m = parameters.mass_total / static_cast<float>(N_edge);
    
    //Jingyi's work 2022-3-3
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            if (cloth.contact_info[index].connecting_spheres.size() > 0)
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

void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint)
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

void rotate_box(bounding_box& box, cgp::inputs_interaction_parameters& inputs)
{
    vec2 const& p1 = inputs.mouse.position.current;
    vec2 const& p0 = inputs.mouse.position.previous;

    bool const event_valid = !inputs.mouse.on_gui;
    bool const click_left = inputs.mouse.click.left;
    bool const click_right = inputs.mouse.click.right;

    if (event_valid) { // If the mouse cursor is not on the ImGui area

        if (click_right)     // Rotation of the camera around its center
            box.manipulator_rotate_trackball(p0, p1);
        else if (click_left) // Translate/Pan the camera in the viewspace plane
            box.manipulator_translate_in_plane(p1 - p0);
    }
}
