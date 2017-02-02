
#define NFACES 4

typedef struct {
  double* vertices_x; // The positions of vertices along x dimension
  double* vertices_y; // The positions of vertices along y dimension
  int* edges_vertices1;        // A vertex adjoining the edge
  int* edges_vertices2;        // The other vertex adjoining the edge

  int* nfaces;        // The number of faces for the cell
  double* density;    // The cell centered density 
  double* energy;     // The cell centered energy

  int* edges;          // nfaces edges for each of the cells

  int* cells_indirection1;  // One of the two cells adjoined by the edge
  int* cells_indirection2;  // The other cell adjoined by the edge

} UnstructuredMesh;

// Build an unstructured mesh
void build_unstructured_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int local_nx, const int local_ny);

