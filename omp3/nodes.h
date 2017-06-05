
// Initialises the CG solver
double initialise_cg(
    const int nx, const int ny, const double dt, const double conductivity,
    const double heat_capacity, double* p, double* r, 
    const double* temperature, const double* volume, const double* b,
    const double* rho, const int* cells_edges, const int* edge_vertex0,
    const int* edge_vertex1, const double* vertices_x, const double* vertices_y,
    const double* cell_centroids_x, const double* cell_centroids_y, 
    const int* edges_cells);

// Calculate the RHS including the unstructured correction term
void calculate_rhs(
    const int nx, const int ny, const double heat_capacity,
    const double conductivity, const double dt, const double* volume, 
    const double* rho, const double* temperature, double* b, 
    const int* edge_vertex0, const int* edge_vertex1,
    const double* cell_centroids_x, const double* cell_centroids_y,
    const double* vertices_x, const double* vertices_y,
    const int* cells_edges, const int* edges_cells);

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, double* p, double* Ap, const double dt, 
    const double conductivity, const double heat_capacity, 
    const double* temperature, const double* volume, const double* rho, 
    const int* cells_edges, const int* edge_vertex0, const int* edge_vertex1, 
    const double* vertices_x, const double* vertices_y,
    const double* cell_centroids_x, const double* cell_centroids_y, 
    const int* edges_cells);

// Updates the current guess using the calculated alpha
double calculate_new_r2(
    int nx, int ny, double alpha, double* x, double* p, double* r, double* Ap);

// Updates the conjugate from the calculated beta and residual
void update_conjugate(
    const int nx, const int ny, const double beta, const double* r, double* p);

// Prints the vector to std out
void print_vec(
    const int nx, const int ny, double* a);

