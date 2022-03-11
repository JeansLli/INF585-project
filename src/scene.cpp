#include "scene.hpp"


using namespace cgp;


void scene_structure::display()
{
	// Basics common elements
	// ***************************************** //
	timer.update();

	remove_sphere();

	environment.light = environment.camera.position();
	if (gui.display_frame)
		draw(global_frame, environment);
	
	
	// Elements of the scene: Obstacles (floor, sphere), and fixed position
	// ***************************************** //
	//draw(obstacle_floor, environment);
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
		simulation_collision_detection(cloth, spheres, box);
		simulation_update_spheres(cloth, spheres, parameters.dt / N_step);
		update_cloth_constraints(cloth, spheres);
		simulation_compute_force(cloth, parameters);

		// One step of numerical integration
		simulation_numerical_integration(cloth, parameters, parameters.dt/N_step);

		// Apply the positional (and velocity) constraints
		simulation_apply_constraints(cloth, constraint);

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

	// Display spheres
	for (size_t k = 0; k < spheres.size(); ++k)
	{
		sphere_structure const& sphere = spheres[k];
		sphere_drawable.transform.translation = sphere.p;
		sphere_drawable.transform.scaling = sphere.r;
		sphere_drawable.shading.color = sphere.c;
		draw(sphere_drawable, environment);
	}

	// Display the box in which the particles should stay
	cgp::buffer<vec3> cube_wireframe_data_new;
	for (auto p : cube_wireframe_data) {
		cube_wireframe_data_new.push_back(box.orientation * box.scale * p + box.center);
	}
	cube_wireframe.update(cube_wireframe_data_new);
	draw(cube_wireframe, environment);
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

	//obstacle_floor.initialize(mesh_primitive_quadrangle({ -1.5f,-1.5f,0 }, { -1.5f,1.5f,0 }, { 1.5f,1.5f,0 }, { 1.5f,-1.5f,0 }));
	//obstacle_floor.texture = opengl_load_texture_image("assets/wood.jpg");
	//obstacle_floor.transform.translation = { 0,0,constraint.ground_z };

	sphere_drawable.initialize(mesh_primitive_sphere());
	
	sphere_fixed_position.initialize(mesh_primitive_sphere());
	sphere_fixed_position.transform.scaling = 0.02f;
	sphere_fixed_position.shading.color = { 0,0,1 };

	cloth_texture = opengl_load_texture_image("assets/cloth.jpg");
	initialize_cloth(gui.N_sample_edge);
	initialize_spheres();
	initialize_box();

	cube_wireframe.initialize(cube_wireframe_data, "cube wireframe");
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

	if (update) {
		cloth.precompute_neighbor(parameters.resolution1, parameters.resolution2);
	}

	ImGui::SliderFloat("Wind direction.x", &parameters.wind.direction.x, -1.0f, 1.0f);
	ImGui::SliderFloat("Wind direction.y", &parameters.wind.direction.y, -1.0f, 1.0f);
	ImGui::SliderFloat("Wind direction.z", &parameters.wind.direction.z, -1.0f, 1.0f);

	bool update_sphere_pos = false;
	update_sphere_pos |= ImGui::SliderFloat("Sphere pos.x", &parameters.sphere_pos.x, -5.0f, 5.0f);
	update_sphere_pos |= ImGui::SliderFloat("Sphere pos.y", &parameters.sphere_pos.y, -5.0f, 5.0f);
	update_sphere_pos |= ImGui::SliderFloat("Sphere pos.z", &parameters.sphere_pos.z, 0.0f, 10.0f);
	if (update_sphere_pos) {
		spheres[cur_sphere_index].p = parameters.sphere_pos;
	}

	ImGui::SliderFloat("Emit velocity.x", &parameters.v_emit.x, 0.0f, 10.0f);
	ImGui::SliderFloat("Emit velocity.y", &parameters.v_emit.y, 0.0f, 10.0f);
	ImGui::SliderFloat("Emit velocity.z", &parameters.v_emit.z, -10.0f, 0.0f);
	
	ImGui::SliderFloat("Box scale", &box.scale, 0.1f, 1.0f);
	ImGui::SliderFloat("Box pos.x", &box.center.x, -10.0f, 10.0f);
	ImGui::SliderFloat("Box pos.y", &box.center.y, -10.0f, 10.0f);
	ImGui::SliderFloat("Box pos.z", &box.center.z, -10.0f, 10.0f);

	ImGui::Spacing(); ImGui::Spacing();
	
	reset |= ImGui::SliderInt("Cloth samples", &gui.N_sample_edge, 4, 80);

	ImGui::Spacing(); ImGui::Spacing();
	reset |= ImGui::Button("Restart");
	if (reset) {
		initialize_spheres();
		initialize_box();
		initialize_cloth(gui.N_sample_edge);
		simulation_running = true;
	}
}

void scene_structure::initialize_box()
{
	// initialize bounding box
	box.initialize(environment.camera);
	box.center = { 0.0f, 1.54f, 0.515f }; // magic number :), try it one by one
	box.scale = 0.35f; // magic number :), try it one by one
}

void scene_structure::initialize_spheres()
{
	spheres.clear();
	parameters.v_emit = { 0.0f, 1.34f, -4.897f }; // magic number :), of course I spend time on trying it
	parameters.sphere_pos = { 0.0f, 0.0f, 1.5f };
	create_sphere();
}

void scene_structure::create_sphere()
{
	sphere_structure sphere;
	sphere.p = parameters.sphere_pos;
	sphere.v = { 0.0f, 0.0f, 0.0f };
	sphere.c = { 1.0f, 1.0f, 1.0f }; // color
	sphere.r = 0.05f;
	sphere.m = 1;
	sphere.can_move = false;
	cur_sphere_index = spheres.size();
	spheres.push_back(sphere);
}

void scene_structure::emit_sphere()
{
	sphere_structure& sphere = spheres[cur_sphere_index];
	sphere.v = parameters.v_emit;
	sphere.can_move = true;
}

void scene_structure::remove_sphere()
{
	if (spheres.size() > 0) {
		auto remove_it = std::remove_if(spheres.begin(), spheres.end(), [&](sphere_structure e) { return e.p.z <= -3.0f; });
		for (auto it = remove_it; it < spheres.end(); it++) {
			for (auto index : it->connecting_particles) {
				size_t sphere_index = it - spheres.begin();
				cloth_contact_info& info = cloth.contact_info[index];
				// remove sphere index
				auto remove_it2 = std::remove_if(info.connecting_spheres.begin(), info.connecting_spheres.end(), [&](size_t e) { return e == sphere_index; });
				for (auto it2 = remove_it2; it2 < info.connecting_spheres.end(); it2++) {
					info.connecting_spheres.erase(it2);
				}
			}
			spheres.erase(it);
		}
	}
}

void scene_structure::keyboard_callback(const cgp::inputs_interaction_parameters& inputs) {
	if (inputs.keyboard.space) {
		emit_sphere();
		create_sphere();
	}
}
