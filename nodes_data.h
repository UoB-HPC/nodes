#ifndef __NODESDATAHDR
#define __NODESDATAHDR

#include "../mesh.h"

#define NEDGES 4
#define NNEIGHBOURS 6   // This is max size required - for 3d
#define ARCH_ROOT_PARAMS "../arch.params"
#define NODES_PARAMS "nodes.params"
#define EPS 1.0e-10
#define NNEIGHBOURS_STENCIL 5
#define NCELLS_PER_EDGE 2

enum { NORTH_STENCIL, EAST_STENCIL, SOUTH_STENCIL, WEST_STENCIL };
enum { BOTTOM, LEFT, RIGHT, TOP };

typedef struct {

  double heat_capacity;
  double conductivity;

  double* b;

} NodesData;

typedef struct {

  // Handles unstructured mesh
  double* vertices_x;
  double* vertices_y;
  double* cell_centroids_x;
  double* cell_centroids_y;
  double* volume;
  int* cells_vertices;
  int* edge_vertex0;
  int* edge_vertex1;
  int* cells_edges;
  int* edges_cells;

  int nedges;

} UnstructuredMesh;

// Initialises the data specific to the nodes application
void initialise_nodes_data(
    const int nx, const int ny, NodesData* nodes_data, const char* nodes_params);

// Build an unstructured mesh
void initialise_unstructured_quad_mesh_2d(
    UnstructuredMesh* unstructured_mesh, Mesh* mesh);

#endif

