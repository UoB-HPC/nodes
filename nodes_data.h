
#define NEDGES 4
#define ARCH_ROOT_PARAMS "../arch.params"
#define NODES_PARAMS "nodes.params"

#define EPS 1.0e-10
#define NNEIGHBOURS_STENCIL 5

enum { NORTH_STENCIL, EAST_STENCIL, SOUTH_STENCIL, WEST_STENCIL };
enum { BOTTOM, LEFT, RIGHT, TOP };

typedef struct {
  int nneighbours;
  int* neighbours_ii;
  int* neighbours_jj;

  double heat_capacity;
  double conductivity;

} NodesData;

typedef struct {

  double* vertices_x;
  double* vertices_y;
  double* cell_centroids_x;
  double* cell_centroids_y;
  int* cells_vertices;
  int* edge_vertex0;
  int* edge_vertex1;
  int* cells_edges;

  double width;
  double height;
  int nedges;

} UnstructuredMesh;

void initialise_nodes_data(
    NodesData* nodes_data, const int global_nx, const int global_ny, 
    const int nx, const int ny, const int nneighbours, const char* nodes_params);

