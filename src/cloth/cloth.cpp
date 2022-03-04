#include "cloth.hpp"

using namespace cgp;


void cloth_structure::initialize(int N_samples_edge_arg)
{
    assert_cgp(N_samples_edge_arg > 3, "N_samples_edge=" + str(N_samples_edge_arg) + " should be > 3");

    position.clear();
    normal.clear();
    velocity.clear();
    force.clear();
    

    position.resize(N_samples_edge_arg, N_samples_edge_arg);
    normal.resize(N_samples_edge_arg, N_samples_edge_arg);
    velocity.resize(N_samples_edge_arg, N_samples_edge_arg);
    force.resize(N_samples_edge_arg, N_samples_edge_arg);
    

    float const z0 = 1.0f;
    mesh const cloth_mesh = mesh_primitive_grid({ -0.5f,0,z0 }, { 0.5f,0,z0 }, { 0.5f,1,z0 }, { -0.5f,1,z0 }, N_samples_edge_arg, N_samples_edge_arg).fill_empty_field();
    position = grid_2D<vec3>::from_buffer(cloth_mesh.position, N_samples_edge_arg, N_samples_edge_arg);
    normal = grid_2D<vec3>::from_buffer(cloth_mesh.normal, N_samples_edge_arg, N_samples_edge_arg);
    triangle_connectivity = cloth_mesh.connectivity;

    neighbor_cross_map.resize(position.size());
    neighbor_diagonal_map.resize(position.size());

    //Jingyi's work 2022-3-4
    initialize_contact_sphere();
    //Jingyi's work done

    precompute_neighbor(2, 1);
}

void cloth_structure::update_normal()
{
    normal_per_vertex(position.data, triangle_connectivity, normal.data);
}

int cloth_structure::N_samples_edge() const
{
    return position.dimension.x;
}

void cloth_structure::initialize_contact_sphere(){
    contact_info.resize(position.size());
    for (int i = 0; i < position.size(); ++i) {
        contact_info[i].is_contact = false;
    }
}

void cloth_structure::initialize_surface_normal(cgp::vec3 normal) {
    surface_normal = normal;
}

void cloth_structure::precompute_neighbor(int resolution1, int resolution2) 
{
    size_t const N_edge = N_samples_edge(); // number of vertices in one dimension of the grid
    for (int ku = 0; ku < N_edge; ++ku) {
        for (int kv = 0; kv < N_edge; ++kv) {
            int index = position.index_to_offset(ku, kv);
            neighbor_cross_map[index].clear();
            neighbor_diagonal_map[index].clear();
            // neighbor and bending springs
            for (int offset = 1; offset <= resolution1; offset++) {
                if ((ku - offset) >= 0) {
                    // left
                    neighbor_cross_map[index].push_back(int3{ ku - offset, kv, offset });
                }
                if ((ku + offset) < N_edge) {
                    // right
                    neighbor_cross_map[index].push_back(int3{ ku + offset, kv, offset });
                }
                if ((kv + offset) < N_edge) {
                    // top
                    neighbor_cross_map[index].push_back(int3{ ku, kv + offset, offset });
                }
                if ((kv - offset) >= 0) {
                    // bottom
                    neighbor_cross_map[index].push_back(int3{ ku, kv - offset, offset });
                }
            }
            // shearing springs
            for (int offset = 1; offset <= resolution2; offset++) {
                if ((ku - offset >= 0) && (kv + offset < N_edge)) {
                    // left-top
                    neighbor_diagonal_map[index].push_back(int3{ ku - offset, kv + offset, offset });
                }
                if ((ku - offset >= 0) && (kv - offset >= 0)) {
                    // left-bottom
                    neighbor_diagonal_map[index].push_back(int3{ ku - offset, kv - offset, offset });
                }
                if ((ku + offset) < N_edge && (kv + offset) < N_edge) {
                    // right-top
                    neighbor_diagonal_map[index].push_back(int3{ ku + offset, kv + offset, offset });
                }
                if ((ku + offset) < N_edge && (kv - offset) >= 0) {
                    // right-bottom
                    neighbor_diagonal_map[index].push_back(int3{ ku + offset, kv - offset, offset });
                }
            }
        }
    }
}

void cloth_structure_drawable::initialize(int N_samples_edge)
{
    mesh const cloth_mesh = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N_samples_edge, N_samples_edge).fill_empty_field();

    drawable.clear();
    drawable.initialize(cloth_mesh, "Cloth");
    drawable.shading.phong.specular = 0.0f;
}

void cloth_structure_drawable::update(cloth_structure const& cloth)
{
    drawable.update_position(cloth.position.data);
    drawable.update_normal(cloth.normal.data);
}
