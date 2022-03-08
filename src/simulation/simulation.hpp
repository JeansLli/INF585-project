#pragma once

#include "cgp/cgp.hpp"
#include "../cloth/cloth.hpp"
#include "../constraint/constraint.hpp"


struct bounding_box
{
    cgp::buffer<cgp::vec3> fns; // face normal
    cgp::buffer<cgp::vec3> fns_axis; // face normal is look at, then this value store its right and up
    cgp::buffer<cgp::vec3> ps; // points on faces

    float scale;
    float width; // original width for each edge of cube
    cgp::rotation_transform orientation;
    cgp::vec3 center;

    cgp::camera_around_center* camera;

    void initialize(cgp::camera_around_center& cam)
    {
        camera = &cam;

        orientation = cgp::rotation_transform();
        center = cgp::vec3(0, 0, 0);

        width = 2;

        // bounding
        fns.push_back(cgp::vec3(0, 1, 0));
        ps.push_back(cgp::vec3(0, -1, 0));
        fns_axis.push_back(cgp::vec3(1, 0, 0)); // right
        fns_axis.push_back(cgp::vec3(0, 0, 1)); // up

        fns.push_back(cgp::vec3(0, -1, 0));
        ps.push_back(cgp::vec3(0, 1, 0));
        fns_axis.push_back(cgp::vec3(-1, 0, 0)); // right
        fns_axis.push_back(cgp::vec3(0, 0, 1)); // up

        fns.push_back(cgp::vec3(1, 0, 0));
        ps.push_back(cgp::vec3(-1, 0, 0));
        fns_axis.push_back(cgp::vec3(0, -1, 0));
        fns_axis.push_back(cgp::vec3(0, 0, 1));

        fns.push_back(cgp::vec3(-1, 0, 0));
        ps.push_back(cgp::vec3(1, 0, 0));
        fns_axis.push_back(cgp::vec3(0, 1, 0));
        fns_axis.push_back(cgp::vec3(0, 0, 1));

        fns.push_back(cgp::vec3(0, 0, 1)); // bottom
        ps.push_back(cgp::vec3(0, 0, -1));
        fns_axis.push_back(cgp::vec3(0, 1, 0));
        fns_axis.push_back(cgp::vec3(1, 0, 0));

        //fns.push_back(cgp::vec3(0, 0, -1)); // top face
        //ps.push_back(cgp::vec3(0, 0, 1)); // point on top face
        //fns_axis.push_back(cgp::vec3(0, -1, 0));
        //fns_axis.push_back(cgp::vec3(1, 0, 0));
    }

    void manipulator_rotate_trackball(cgp::vec2 const& p0, cgp::vec2 const& p1)
    {
        /*cgp::rotation_transform const r = trackball_rotation(p0, p1);
        orientation = orientation * inverse(r);*/

        cgp::vec2 t2d = p1 - p0;
        const float d = norm(t2d);
        if (d > 1e-4f) {
            float theta = t2d.x;
            float phi = t2d.y;
            cgp::rotation_transform r1 = cgp::rotation_transform::from_axis_angle(camera->up(), theta);
            cgp::rotation_transform r2 = cgp::rotation_transform::from_axis_angle(-camera->right(), phi);
            orientation *= (r2 * r1);
        }
    }

    void manipulator_translate_in_plane(cgp::vec2 const& tr)
    {
        center += translation_in_plane(tr, camera->orientation());
    }

    cgp::vec3 get_surface_normal(int i) const
    {
        cgp::vec3 n = fns[i];
        return normalize(orientation * n);
    }

    cgp::vec3 get_point_on_surface(int i) const
    {
        cgp::vec3 p = ps[i];
        return orientation * scale * p + center;
    }

    cgp::vec3 get_right_vector(int i) const
    {
        cgp::vec3 v = fns_axis[i * 2];
        return normalize(orientation * v);
    }

    cgp::vec3 get_up_vector(int i) const
    {
        cgp::vec3 v = fns_axis[i * 2 + 1];
        return normalize(orientation * v);
    }

    float get_width() const
    {
        return scale * width;
    }
};

void rotate_box(bounding_box& box, cgp::inputs_interaction_parameters& inputs);

struct particle_structure
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed

    cgp::vec3 c; // Color
    float r;     // Radius
    float m;     // mass
};

struct sphere_structure :public particle_structure {
    bool can_move;
    bool detect_cloth_collision; // whether it needs detection of collision between sphere and cloth
    cgp::vec3 p_at_contact; // p while first contact
    std::vector<size_t> connecting_particles; // connecting particles of cloth
};

struct simulation_parameters
{
    float dt = 0.02f;        // time step for the numerical integration
    float mass_total = 1.0f; // total mass of the cloth
    float K = 10.0f;         // stiffness parameter
    float mu = 40.0f;        // damping parameter
    int resolution1 = 2;     // neighbor&bending resolution
    int resolution2 = 1;     // shearing resolution
    cgp::vec3 sphere_pos;   // start pos
    cgp::vec3 v_emit;  // velocity of emit

    //  Wind magnitude and direction
    struct {
        float magnitude = 0.0f;
        cgp::vec3 direction = { 0,-1,0 };
    } wind;
};


// Fill the forces in the cloth given the position and velocity
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters);

// Perform 1 step of a semi-implicit integration with time step dt
void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt);

// Apply the constraints (fixed position, obstacles) on the cloth position and velocity
void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint);

// Helper function that tries to detect if the simulation diverged 
bool simulation_detect_divergence(cloth_structure const& cloth);

void simulation_collision_detection(cloth_structure& cloth, std::vector<sphere_structure>& spheres, bounding_box const& box);

void collision_sphere_sphere(std::vector<sphere_structure>& spheres);

void collision_sphere_cloth(cloth_structure& cloth, std::vector<sphere_structure>& spheres);

void simulation_update_spheres(cloth_structure& cloth, std::vector<sphere_structure>& spheres, float dt);

void update_cloth_constraints(cloth_structure& cloth, std::vector<sphere_structure>& spheres);


