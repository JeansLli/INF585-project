#include "scene.hpp"


using namespace cgp;


void scene_structure::display()
{
	// Basics common elements
	// ***************************************** //
	timer.update();
	environment.light = environment.camera.position();
	if (gui.display_frame)
		draw(global_frame, environment);

	
	
	// Elements of the scene: Obstacles (floor, sphere), and fixed position
	// ***************************************** //
	draw(obstacle_floor, environment);
	draw(obstacle_sphere, environment);
	for (auto const& c : constraint.fixed_sample)
	{
		sphere_fixed_position.transform.translation = c.second.position;
		draw(sphere_fixed_position, environment);
	}


	// Simulation of the cloth
	// ***************************************** //
	int const N_step = 1; // Adapt here the number of intermediate simulation steps (ex. 5 intermediate steps per frame)
	for (int k_step = 0; simulation_running ==true && k_step < N_step; ++k_step)
	{
		// Update the forces on each particle
		simulation_compute_force(cloth, parameters);

		// One step of numerical integration
		simulation_numerical_integration(cloth, falling_sphere, parameters, parameters.dt/N_step);

		// Apply the positional (and velocity) constraints
		simulation_apply_constraints(cloth, falling_sphere, constraint, parameters);


		// Check if the simulation has not diverged - otherwise stop it
		bool const simulation_diverged = simulation_detect_divergence(cloth);
		if (simulation_diverged) {
			std::cout << "\n *** Simulation has diverged ***" << std::endl;
			std::cout << " > The simulation is stoped" << std::endl;
			simulation_running = false;
		}
	}


	// Cloth display
	// ***************************************** //
	
	// Prepare to display the updated cloth
	cloth.update_normal();        // compute the new normals
	cloth_drawable.update(cloth); // update the positions on the GPU

	// Display the cloth
	draw(cloth_drawable, environment);
	if(gui.display_wireframe)
		draw_wireframe(cloth_drawable, environment);

	// Display falling sphere
	falling_sphere_drawable.transform.translation = falling_sphere.p;
	falling_sphere_drawable.transform.scaling = falling_sphere.r;
	falling_sphere_drawable.shading.color = falling_sphere.c;
	draw(falling_sphere_drawable, environment);
}


// Compute a new cloth in its initial position (can be called multiple times)
void scene_structure::initialize_cloth(int N_sample)
{
	cloth.initialize(N_sample);
	cloth_drawable.initialize(N_sample);
	cloth_drawable.drawable.texture = cloth_texture;

	constraint.fixed_sample.clear();
	constraint.add_fixed_position(0, 0, cloth);
	constraint.add_fixed_position(0, N_sample - 1, cloth);
	constraint.add_fixed_position(N_sample - 1, 0, cloth);
	constraint.add_fixed_position(N_sample - 1, N_sample - 1, cloth);
}


void scene_structure::initialize()
{
	global_frame.initialize(mesh_primitive_frame(), "Frame");
	environment.camera.look_at({ 3.0f,2.0f,2.0f }, { 0,0,0 }, { 0,0,1 });

	obstacle_floor.initialize(mesh_primitive_quadrangle({ -1.5f,-1.5f,0 }, { -1.5f,1.5f,0 }, { 1.5f,1.5f,0 }, { 1.5f,-1.5f,0 }));
	obstacle_floor.texture = opengl_load_texture_image("assets/wood.jpg");
	obstacle_floor.transform.translation = { 0,0,constraint.ground_z };

	obstacle_sphere.initialize(mesh_primitive_sphere());
	obstacle_sphere.transform.translation = constraint.sphere.center;
	obstacle_sphere.transform.scaling = constraint.sphere.radius;
	obstacle_sphere.shading.color = { 1,0,0 };

	falling_sphere_drawable.initialize(mesh_primitive_sphere());
	
	sphere_fixed_position.initialize(mesh_primitive_sphere());
	sphere_fixed_position.transform.scaling = 0.02f;
	sphere_fixed_position.shading.color = { 0,0,1 };


	cloth_texture = opengl_load_texture_image("assets/cloth.jpg");
	initialize_cloth(gui.N_sample_edge);

	initialize_spheres();
}

void scene_structure::display_gui()
{
	bool reset = false;

	ImGui::Text("Display");
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
	ImGui::Checkbox("Texture Cloth", &cloth_drawable.drawable.shading.use_texture);

	ImGui::Spacing(); ImGui::Spacing();

	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Time step", &parameters.dt, 0.01f, 0.03f);
	ImGui::SliderFloat("Stiffness", &parameters.K, 1.0f, 80.0f, "%.3f", 2.0f);
	ImGui::SliderFloat("Wind magnitude", &parameters.wind.magnitude, 0, 60, "%.3f", 2.0f);
	ImGui::SliderFloat("Damping", &parameters.mu, 1.0f, 100.0f);
	ImGui::SliderFloat("Mass", &parameters.mass_total, 0.2f, 5.0f, "%.3f", 2.0f);

	bool update = false;
	update |= ImGui::SliderInt("Neighbor&Bending Resolution", &parameters.resolution1, 2, 10);
	update |= ImGui::SliderInt("Shearing Sprins Resolution", &parameters.resolution2, 1, 10);

	ImGui::SliderFloat("Wind direction.x", &parameters.wind.direction.x, -1.0f, 1.0f);
	ImGui::SliderFloat("Wind direction.y", &parameters.wind.direction.y, -1.0f, 1.0f);
	ImGui::SliderFloat("Wind direction.z", &parameters.wind.direction.z, -1.0f, 1.0f);

	if (update) {
		cloth.precompute_neighbor(parameters.resolution1, parameters.resolution2);
	}

	ImGui::Spacing(); ImGui::Spacing();
	
	reset |= ImGui::SliderInt("Cloth samples", &gui.N_sample_edge, 4, 80);

	ImGui::Spacing(); ImGui::Spacing();
	reset |= ImGui::Button("Restart");
	if (reset) {
		initialize_spheres();
		initialize_cloth(gui.N_sample_edge);
		simulation_running = true;
	}

}

void scene_structure::initialize_spheres()
{
	falling_sphere.p = { 0.0f, 0.5f, 1.5f };
	falling_sphere.v = { 0.0f, 0.0f, 0.0f };
	falling_sphere.c = { 1.0f, 1.0f, 1.0f }; // color
	falling_sphere.r = 0.05f;
	falling_sphere.m = 1;
}
