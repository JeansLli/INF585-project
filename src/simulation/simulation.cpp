#include "simulation.hpp"
#include <math.h>

using namespace cgp;

const vec3 g = { 0,0,-9.81f };
//const float alpha = 0.9;// attenuation for collision
//const float beta = 0.9;// attenuation for collision
const float epsilon = 0.01f;
float t = 0;//time for extension
bool EXTENSION_CLOTH = false;
bool CONTRACTION_CLOTH = false;
float sigma = 3;// the smaller sigma is, the slower the sphere's speed decrease


static vec3 spring_force(const vec3& p_i, const vec3& p_j, float L0, float K)
{
    vec3 const p = p_i - p_j;
    float const L = norm(p);
    vec3 const u = p / L;
    vec3 const F = -K * (L - L0) * u;
    return F;
}


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
            v1 = v1 / 1.2f;
            v2 = v2 / 1.2f;
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


void simulation_collision_detection(cloth_structure &cloth,particle_structure& falling_sphere){
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    const float radius_offset = 0.01f; // because we just check the vertices not triangles, that will lead to vertices on sphere's surface but the middle point get inside into sphere
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            vec3& p = cloth.position(ku, kv);
            vec3 cf = p - falling_sphere.p;
            float cf_len = norm(cf);
            if(cf_len <= (falling_sphere.r + radius_offset)){
                cloth.contact_sphere[index]= true;
            }
            else{
                cloth.contact_sphere[index]=false;
            }
        }
    }         
}

void fall_sphere_update(particle_structure& falling_sphere, float dt){
    if(EXTENSION_CLOTH==false && CONTRACTION_CLOTH==false){
        //std::cout << "fall freely"<<std::endl;

        vec3 const f = 0.05f * g;

        falling_sphere.v = (1 - 0.9f * dt) * falling_sphere.v + dt * f;
        falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
        //std::cout << "no extension, v=" << falling_sphere.v << std::endl;

    }
    else if (EXTENSION_CLOTH==true && CONTRACTION_CLOTH==false){
        //std::cout << "extension"<<std::endl;
        if(falling_sphere.BALL_STOP){
            return;
        }
        if(norm(falling_sphere.v)<1e-7){
            t = t - dt;
            falling_sphere.v = -falling_sphere.v;
            falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
            EXTENSION_CLOTH = false;
            CONTRACTION_CLOTH = true;
        }else{
            float gaussion_decrease_last = exp(-((t-dt)*sigma)*((t-dt)*sigma)/2);
            float gaussion_decrease = exp(-(t*sigma)*(t*sigma)/2);
            //std::cout << "gaussion_decrease="<<gaussion_decrease << std::endl;
            falling_sphere.v = falling_sphere.v * gaussion_decrease / gaussion_decrease_last; // v_k = v_k-1 * g(t_k) / g_(t_k-1)
            //std::cout << "extension, v=" << falling_sphere.v << std::endl;
            falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
            //std::cout << "extension, p=" << falling_sphere.p << std::endl;
            //std::cout << "extension, move p=" << dt * falling_sphere.v << std::endl;
            t = t + dt;
        }

    }
    else if(EXTENSION_CLOTH==false && CONTRACTION_CLOTH==true){
        //std::cout << "contrction"<<std::endl;
        float gaussion_decrease_last = exp(-((t+dt)*sigma)*((t+dt)*sigma)/2);
        float gaussion_decrease = exp(-(t*sigma)*(t*sigma)/2);
        //std::cout << "gaussion_decrease="<<gaussion_decrease << std::endl;
        falling_sphere.v = falling_sphere.v * gaussion_decrease / gaussion_decrease_last; // v_k = v_k-1 * g(t_k) / g_(t_k-1)
        //std::cout << "extension, v=" << falling_sphere.v << std::endl;
        falling_sphere.p = falling_sphere.p + dt * falling_sphere.v;
        t = t - dt;
        if(t<1e-9){
            EXTENSION_CLOTH = false;
            CONTRACTION_CLOTH = false;
        }
    }
}

void update_cloth_constraints(cloth_structure& cloth, particle_structure& falling_sphere){
    size_t const N_edge = cloth.N_samples_edge(); // number of vertices in one dimension of the grid
    const float radius_offset = 0.01f; // because we just check the vertices not triangles, that will lead to vertices on sphere's surface but the middle point get inside into sphere
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            if(cloth.contact_sphere[index]){
                vec3& p = cloth.position(ku, kv);
                
            }
            
            // collide with obstacle spheres
            /*
            vec3 normal = cloth.normal(ku, kv);
            vec3 cp = p - constraint.sphere.center; // center to p
            float cp_len = norm(cp);
            float diff = cp_len - (constraint.sphere.radius + radius_offset);
            if (diff <= 0) {
                //p += cp / cp_len * diff; // move along with center to p direction -> not good, will cause oscillation
                if(dot(normal,cp)>=0){
                    p += normal * std::abs(dot(cp / cp_len * diff, normal)); // move along with center to p direction
                }else{
                    p += -normal * std::abs(dot(cp / cp_len * diff, normal));
                }
            }
            */

            // collide with ground, add epsilon to avoid rendering cloth and plane alternately when their triangle overlap
            //p.z = std::max(constraint.ground_z + epsilon, p.z);

            //vec3 normal = cloth.normal(ku, kv);
            vec3 cf = p - falling_sphere.p;
            float cf_len = norm(cf);
            //vec3 uij = (p - falling_sphere.p) / detection; // direction from pi(cloth particle) to pj(falling_sphere)
            //falling_sphere.r + radius_offset=0.06
            if (cf_len <= (falling_sphere.r + radius_offset)) {
                //Jingyi's work 2022-3-3
                //collision_sphere_sphere(falling_sphere,cloth,ku,kv,radius_offset); //it's for colission

                //It's for Gaussian decay and cloth extension 
                
                //std::cout << "position before="<<p<<std::endl;
                //if(dot(cf,falling_sphere.v)>=0){

                //if(cf.z <=0){ //the sphere's center is above
                    p += (falling_sphere.r + radius_offset - cf_len) * normalize(cf); 
                    if(t<1e-9 && falling_sphere.v.z<0){
                        EXTENSION_CLOTH = true;
                        CONTRACTION_CLOTH = false;
                        //falling_sphere.BALL_STOP = true;
                    }
                    

                    //std::cout << "move distance="<<(falling_sphere.r + radius_offset - cf_len)<<std::endl;
                    //std::cout << "move direction="<<normalize(cf)<<std::endl;
                //}
                //else{
                    //p -= (falling_sphere.r + radius_offset + cf_len) * normalize(cf);
                    //std::cout << "move distance="<<(falling_sphere.r + radius_offset + cf_len)<<std::endl;
                    //std::cout << "move direction="<<-normalize(cf)<<std::endl;
                //}
                //std::cout << "position after="<<p<<std::endl;
                //std::cout << "falling_sphere position"<<falling_sphere.p<<std::endl;
                //std::cout << "falling_sphere radius"<<falling_sphere.r<<std::endl;
                //EXTENSION_CLOTH = true;
                //Jingyi's work done
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
    
    // Springs
    const int N_neighbors = 2;

    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            vec3& f = force(ku, kv);

            // Loops over neighbors
            for (int du = -N_neighbors; du <= N_neighbors; ++du) {
                for (int dv = -N_neighbors; dv <= N_neighbors; ++dv) {
                    if (du != 0 || dv != 0)
                    {
                        // Neighbors indices
                        int const ku_n = ku + du;
                        int const kv_n = kv + dv;
                        float const alpha = float(std::sqrt(du * du + dv * dv)); // rest-length

                        if (ku_n >= 0 && ku_n < N_edge && kv_n >= 0 && kv_n < N_edge)
                            f += spring_force(position(ku, kv), position(ku_n, kv_n), alpha * L0, K / alpha);
                    }
                }
            }

        }

    }
}

void simulation_numerical_integration(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters const& parameters, float dt)
{
    int const N_edge = cloth.N_samples_edge();
    float const m = parameters.mass_total / static_cast<float>(N_edge);
    

    //Jingyi's work 2022-3-3
    //std::cout << "new loopppppppppppp" << std::endl;

    
    
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = cloth.position.index_to_offset(ku, kv);
            vec3& v = cloth.velocity(ku, kv);
            vec3& p = cloth.position(ku, kv);
            if(cloth.contact_sphere[index]){
                std::cout << "contacted particle=" << ku << ","<< kv << std::endl;
                p = p + dt * falling_sphere.v;
            }else{
                vec3 const& f = cloth.force(ku, kv);
                // Standard semi-implicit numerical integration
                v = v + dt * f / m;
                p = p + dt * v;
            }
        }
    }

    //std::cout<<"falling sphere position "<<falling_sphere.p<<std::endl;
    //std::cout<<"falling sphere velocity "<<falling_sphere.v<<std::endl;
    //std::cout<<"t="<<t<<std::endl;


//Jingyi's work done
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