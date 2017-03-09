#pragma once

#include "../mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

  // performs the cg solve, you always want to perform these steps, regardless
  // of the context of the problem etc.
  void solve_unstructured_diffusion_2d(
      const int local_nx, const int local_ny, const int negdes, 
      const double* vertices_x, const double* vertices_y, const int* cells_vertices, 
      const int* nfaces, const double* density, const int* cells_indirection1, 
      const int* cells_indirection2, const int* edges, double* energy);

#ifdef __cplusplus
}
#endif
