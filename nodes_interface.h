#pragma once

#include "../mesh.h"
#include "nodes_data.h"

#ifdef __cplusplus
extern "C" {
#endif

// Solve the unstructured diffusion problem
void solve_unstructured_diffusion_2d(
    const int nx, const int ny, const int pad, Mesh* mesh,
    UnstructuredMesh* unstructured_mesh, const int max_inners, const double dt,
    const double heat_capacity, const double conductivity, double* temperature,
    double* b, double* r, double* p, double* rho, double* Ap, int* end_niters,
    double* end_error, double* reduce_array);

#ifdef __cplusplus
}
#endif
