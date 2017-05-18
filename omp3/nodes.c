#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nodes.h"
#include "../nodes_data.h"
#include "../nodes_interface.h"
#include "../../profiler.h"
#include "../../comms.h"

// Solve the unstructured diffusion problem
void solve_unstructured_diffusion_2d(
    const int nx, const int ny, Mesh* mesh, const int max_inners, const double dt, 
    const double heat_capacity, const double conductivity,
    double* x, double* r, double* p, double* rho, double* s_x, double* s_y, 
    double* Ap, int* end_niters, double* end_error, double* reduce_array,
    const double* edgedx, const double* edgedy, const int nneighbours, 
    int* neighbours_ii, int* neighbours_jj)
{
  // Store initial residual
  double local_old_r2 = initialise_cg(
      nx, ny, dt, heat_capacity, conductivity, p, r, x, rho, s_x, s_y, edgedx, 
      edgedy, nneighbours, neighbours_ii, neighbours_jj);
  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, x, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for(ii = 0; ii < max_inners; ++ii) {

    const double local_pAp = calculate_pAp(nx, ny, s_x, s_y, p, Ap, 
        nneighbours, neighbours_ii, neighbours_jj);
    const double global_pAp = reduce_all_sum(local_pAp);
    const double alpha = global_old_r2/global_pAp;

    const double local_new_r2 = calculate_new_r2(nx, ny, alpha, x, p, r, Ap);
    const double global_new_r2 = reduce_all_sum(local_new_r2);
    const double beta = global_new_r2/global_old_r2;
    handle_boundary_2d(nx, ny, mesh, x, NO_INVERT, PACK);

    // Check if the solution has converged
    if(fabs(global_new_r2) < EPS) {
      global_old_r2 = global_new_r2;
      break;
    }

    update_conjugate(nx, ny, beta, r, p);
    handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);

    // Store the old squared residual
    global_old_r2 = global_new_r2;
  }

  *end_niters = ii;
  *end_error = global_old_r2;
}

// Initialises the CG solver
double initialise_cg(
    const int nx, const int ny, const double dt, const double conductivity,
    const double heat_capacity, double* p, double* r, const double* x, 
    const double* rho, double* s_x, double* s_y, const double* edgedx, 
    const double* edgedy, const int nneighbours, int* neighbours_ii, 
    int* neighbours_jj)
{
  START_PROFILING(&compute_profile);

  // https://inldigitallibrary.inl.gov/sti/3952796.pdf
  // Take the average of the coefficients at the cells surrounding 
  // each face
#pragma omp parallel for
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < (nx+1)-PAD; ++jj) {
      s_x[(ii)*(nx+1)+(jj)] = (dt*conductivity*(rho[(ii)*nx+(jj)]+rho[(ii)*nx+(jj-1)]))/
        (2.0*rho[(ii)*nx+(jj)]*rho[(ii)*nx+(jj-1)]*edgedx[jj]*edgedx[jj]*heat_capacity);
    }
  }
#pragma omp parallel for
  for(int ii = PAD; ii < (ny+1)-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      s_y[(ii)*nx+(jj)] = (dt*conductivity*(rho[(ii)*nx+(jj)]+rho[(ii-1)*nx+(jj)]))/
        (2.0*rho[(ii)*nx+(jj)]*rho[(ii-1)*nx+(jj)]*edgedy[ii]*edgedy[ii]*heat_capacity);
    }
  }

  double initial_r2 = 0.0;
#pragma omp parallel for reduction(+: initial_r2)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int neighbour_index = (ii)*nx*nneighbours+(jj)*nneighbours;
      const int nii = neighbours_ii[neighbour_index+NORTH_STENCIL];
      const int njj = neighbours_jj[neighbour_index+NORTH_STENCIL];
      const int eii = neighbours_ii[neighbour_index+EAST_STENCIL];
      const int ejj = neighbours_jj[neighbour_index+EAST_STENCIL];
      const int sii = neighbours_ii[neighbour_index+SOUTH_STENCIL];
      const int sjj = neighbours_jj[neighbour_index+SOUTH_STENCIL];
      const int wii = neighbours_ii[neighbour_index+WEST_STENCIL];
      const int wjj = neighbours_jj[neighbour_index+WEST_STENCIL];

      r[(ii)*nx+(jj)] = x[(ii)*nx+(jj)] -
        ((s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
          s_x[(eii)*(nx+1)+(ejj)]+s_y[(nii)*nx+(njj)])*x[(ii)*nx+(jj)]
         - s_y[(ii)*nx+(jj)]*x[(sii)*nx+(sjj)]
         - s_x[(ii)*(nx+1)+(jj)]*x[(wii)*nx+(wjj)] 
         - s_x[(eii)*(nx+1)+(ejj)]*x[(eii)*nx+(ejj)]
         - s_y[(nii)*nx+(njj)]*x[(nii)*nx+(njj)]);
      p[(ii)*nx+(jj)] = r[(ii)*nx+(jj)];
      initial_r2 += r[(ii)*nx+(jj)]*r[(ii)*nx+(jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "initialise cg");
  return initial_r2;
}

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, const double* s_x, const double* s_y,
    double* p, double* Ap, const int nneighbours, int* neighbours_ii, 
    int* neighbours_jj)
{
  START_PROFILING(&compute_profile);

  // You don't need to use a matrix as the band matrix is fully predictable
  // from the 5pt stencil
  double pAp = 0.0;
#pragma omp parallel for reduction(+: pAp)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int nii = neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+NORTH_STENCIL];
      const int njj = neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+NORTH_STENCIL];
      const int eii = neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+EAST_STENCIL];
      const int ejj = neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+EAST_STENCIL];
      const int sii = neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+SOUTH_STENCIL];
      const int sjj = neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+SOUTH_STENCIL];
      const int wii = neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+WEST_STENCIL];
      const int wjj = neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+WEST_STENCIL];

      Ap[(ii)*nx+(jj)] = 
        (s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
         s_x[(eii)*(nx+1)+(ejj)]+s_y[(nii)*nx+(njj)])*p[(ii)*nx+(jj)]
        - s_y[(ii)*nx+(jj)]*p[(sii)*nx+(sjj)]
        - s_x[(ii)*(nx+1)+(jj)]*p[(wii)*nx+(wjj)] 
        - s_x[(eii)*(nx+1)+(ejj)]*p[(eii)*nx+(ejj)]
        - s_y[(nii)*nx+(njj)]*p[(nii)*nx+(njj)];
      pAp += p[(ii)*nx+(jj)]*Ap[(ii)*nx+(jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "calculate alpha");
  return pAp;
}

// Updates the current guess using the calculated alpha
double calculate_new_r2(
    int nx, int ny, double alpha, double* x, double* p, double* r, double* Ap)
{
  START_PROFILING(&compute_profile);

  double new_r2 = 0.0;

#pragma omp parallel for reduction(+: new_r2)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      x[(ii)*nx+(jj)] += alpha*p[(ii)*nx+(jj)];
      r[(ii)*nx+(jj)] -= alpha*Ap[(ii)*nx+(jj)];
      new_r2 += r[(ii)*nx+(jj)]*r[(ii)*nx+(jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "calculate new r2");
  return new_r2;
}

// Updates the conjugate from the calculated beta and residual
void update_conjugate(
    const int nx, const int ny, const double beta, const double* r, double* p)
{
  START_PROFILING(&compute_profile);
#pragma omp parallel for
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      p[(ii)*nx+(jj)] = r[(ii)*nx+(jj)] + beta*p[(ii)*nx+(jj)];
    }
  }
  STOP_PROFILING(&compute_profile, "update conjugate");
}

// Prints the vector to std out
void print_vec(
    const int nx, const int ny, double* a)
{
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      printf("%.3e ", a[ii*nx+jj]);
    }
    printf("\n");
  }
}

#if 0
/* At this stage we can consider our unstructured data well initialised.
 *
 * The steps needed to determine the coefficients of the sparse coefficient 
 * matrix are:
 *    (1) Determine the gradients across each face
 *    (2) Determine the second gradients
 *    (3) Calculate the harmonic mean of the density
 *    (4) Initialise the sparse matrix with the data
 * */

/* We are absolutely 100% making the assumption that the grid is not
 * orthogonally partitioned, which means we will use the derivation that
 * incorporates a transformation to the basis of the face normal and 
 * perpendicular vectors, even though in practice the actual vectors
 * will be orthogonal. */

// Determine the gradients across each face
for(int ii = 0; ii < nedges; ++ii) {
  // Fetch the two cells that border the face
  const int cell1 = cells_indirection1[ii];
  const int cell2 = cells_indirection2[ii];

  // TODO: is it correct behaviour here to just leave?
  // Check that we aren't at a boundary, in which case just leave??
  if(cell1 == EDGE || cell2 == EDGE) {
    continue;
  }

#if 0
  grad_edge[] = ;
#endif // if 0
}
#endif // if 0

