#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"
#include "../constraint/constraint.hpp"


struct simulation_parameters
{
    float dt = 0.02f;        // time step for the numerical integration
    float mass_total = 1.0f; // total mass of the cloth
    float K = 10.0f;         // stiffness parameter
    float mu = 40.0f;        // damping parameter
    int resolution1 = 2;     // neighbor&bending resolution
    int resolution2 = 1;     // shearing resolution

    bool need_collision_detect = true; // whether we need collision detection between sphere and cloth
    bool is_connecting = false;  // falling sphere is connecting with cloth
    bool is_extending = false;  // cloth is extending

    //  Wind magnitude and direction
    struct {
        float magnitude = 0.0f;
        cgp::vec3 direction = { 0,-1,0 };
    } wind;

    void reset() {
        need_collision_detect = true;
        is_connecting = false;
        is_extending = false;
    }
};


// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters const& parameters, float dt);

// Apply the constraints (fixed position, obstacles) on the cloth position and velocity
void simulation_apply_constraints(cloth_structure& cloth, particle_structure& falling_sphere, constraint_structure const& constraint, simulation_parameters const& parameters);

// Helper function that tries to detect if the simulation diverged 
bool simulation_detect_divergence(cloth_structure const& cloth);

void simulation_collision_detection(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters& parameters);

//Handle two spheres collision
void collision_sphere_sphere(particle_structure& sphere, cloth_structure& cloth, int ku, int kv, float radius_offset);

void fall_sphere_update(cloth_structure& cloth, particle_structure& falling_sphere, simulation_parameters& parameters, float dt);

void update_cloth_constraints(cloth_structure& cloth, particle_structure& falling_sphere);