
#define NFACES 4
#define ARCH_ROOT_PARAMS "../arch.params"
#define NODES_PARAMS "nodes.params"

#define EPS 1.0e-10
#define NNEIGHBOURS_STENCIL 5

enum { NORTH_STENCIL, EAST_STENCIL, SOUTH_STENCIL, WEST_STENCIL };

typedef struct {
  int nneighbours;
  int* neighbours_ii;
  int* neighbours_jj;

  double heat_capacity;
  double conductivity;

} NodesData;

void initialise_nodes_data(
    NodesData* nodes_data, const int nx, const int ny, const int nneighbours, 
    const char* nodes_params);

#if 0
typedef struct {
  double* values;
  int* cols;
  int* rows;
} csr_matrix;

typedef struct {
  double* vertices_x; // The positions of vertices along x dimension
  double* vertices_y; // The positions of vertices along y dimension

  int ncell_vertices; // The sum of nvertices for every cell
  int* cells_vertices; // Vertices surrounding cell in clockwise order

  int* nfaces;        // The number of faces for the cell
  double* density;    // The cell centered density
  double* energy;     // The cell centered energy

  int nedges;      // The number of edges that exist on the mesh
  int* edges;          // nfaces edges for each of the cells

  int* cells_indirection1;  // One of the two cells adjoined by the edge
  int* cells_indirection2;  // The other cell adjoined by the edge

  double* cell_centers_x; // The x dimension of the center of the cell by index
  double* cell_centers_y; // The y dimension of the center of the cell by index

  double* grad_edge;

  double width;
  double height;

  double conductivity;
  double heat_capacity;
  int max_inner_iterations;

} UnstructuredMesh;

// Build an unstructured mesh
void build_unstructured_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int local_nx, const int local_ny);
#endif // if 0

