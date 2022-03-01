#pragma once

#include "cgp/cgp.hpp"

#include "cloth/cloth.hpp"
#include "simulation/simulation.hpp"

// The element of the GUI that are not already stored in other structures
struct gui_parameters {
	bool display_frame     = false;
	bool display_wireframe = false;
	int N_sample_edge    = 20;  // number of samples of the cloth (the total number of vertices is N_sample_edge^2)
};



// The structure of the custom scene
struct scene_structure {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //

	cgp::mesh_drawable global_frame;          // The standard global frame
	cgp::scene_environment_basic environment; // Standard environment controler
	gui_parameters gui;                       // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	// Display of the obstacles and constraints
	cgp::mesh_drawable obstacle_floor;
	cgp::mesh_drawable obstacle_sphere;
	cgp::mesh_drawable sphere_fixed_position;

	// TODO: create particles to rendering at a time[YUAN test]
	cgp::mesh_drawable falling_sphere_drawable;

	// Cloth related structures
	cloth_structure cloth;                     // The values of the position, velocity, forces, etc, stored as a 2D grid
	cloth_structure_drawable cloth_drawable;   // Helper structure to display the cloth as a mesh
	simulation_parameters parameters;          // Stores the parameters of the simulation (stiffness, mass, damping, time step, etc)
	constraint_structure constraint;           // Handle the parameters of the constraints (fixed vertices, floor and sphere)

	// Helper variables
	bool simulation_running = true;   // Boolean indicating if the simulation should be computed
	GLuint cloth_texture;             // Storage of the texture ID used for the cloth

	//YUAN test, todo, click space to launch sphere
	particle_structure falling_sphere;


	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void display();     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop

	void initialize_cloth(int N_sample); // Recompute the cloth from scratch

	void initialize_spheres();
};





