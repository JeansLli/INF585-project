# INF585-project
This is project is based on code-base from computer animation course, to implement the simulation of cloth, sphere, box and their collision between themselves.

What this project has done:
1. multiple spheres collide with cloth
2. collision between spheres
3. spheres collide with box
4. rotation, translation of box from manipulation
5. wind effect on cloth



Spheres collide with cloth main procedure:

1. simulation_collision_detection(cloth, spheres, box);
2. simulation_update_spheres(cloth, spheres, parameters.dt / N_step);
3. update_cloth_constraints(cloth, spheres);
4. simulation_compute_force(cloth, parameters);
5. simulation_numerical_integration(cloth, parameters, parameters.dt/N_step);
6. simulation_apply_constraints(cloth, constraint);
