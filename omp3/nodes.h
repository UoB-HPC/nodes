
// Initialises the CG solver
double initialise_cg(
    const int nx, const int ny, const double dt, const double conductivity,
    const double heat_capacity, double* p, double* r, const double* x, 
    const double* rho, double* s_x, double* s_y, const double* edgedx, 
    const double* edgedy, const int nneighbours, int* neighbours_ii, 
    int* neighbours_jj);

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, const double* s_x, const double* s_y,
    double* p, double* Ap, const int nneighbours, int* neighbours_ii, int* neighbours_jj);

// Updates the current guess using the calculated alpha
double calculate_new_r2(
    int nx, int ny, double alpha, double* x, double* p, double* r, double* Ap);

// Updates the conjugate from the calculated beta and residual
void update_conjugate(
    const int nx, const int ny, const double beta, const double* r, double* p);

// Prints the vector to std out
void print_vec(
    const int nx, const int ny, double* a);

