#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nodes.h"
#include "../nodes_interface.h"
#include "../../profiler.h"
#include "../../comms.h"

#define ind0 (ii*nx + jj)
#define ind1 (ii*(nx+1) + jj)

/* At the moment we have the shell of the CG solver and actually
 * just going to start with a structured cartesian grid, with the
 * pretense that it is unstructured. */


/* The actual solve will need to be performed by a more general CG
 * implementation that utilises sparse matrix storage i.e. CSR of some form.
 *
 * As the actual matrix will no longer be a predictable band matrix by
 * any general formulation, although in our case it will from the 
 * perspective that the grid will be cartesian in spite of it's unstructured
 * nature. */

// performs the cg solve, you always want to perform these steps, regardless
// of the context of the problem etc.
void solve_unstructured_diffusion_2d(
    const int local_nx, const int local_ny, const int nedges, 
    const double* vertices_x, const double* vertices_y, const int* cells_vertices, 
    const int* nfaces, const double* density, const int* cells_indirection1, 
    const int* cells_indirection2, const int* edges, double* energy)
{
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
  
#if 0
  // Store initial residual
  double local_old_r2 = initialise_cg(
      nx, ny, dt, p, r, x, rho, s_x, s_y, edgedx, edgedy);
  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, x, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for(ii = 0; ii < MAX_INNER_ITERATIONS; ++ii) {

    const double local_pAp = calculate_pAp(nx, ny, s_x, s_y, p, Ap);
    const double global_pAp = reduce_all_sum(local_pAp);
    const double alpha = global_old_r2/global_pAp;

    const double local_new_r2 = calculate_new_r2(nx, ny, alpha, x, p, r, Ap);
    const double global_new_r2 = reduce_all_sum(local_new_r2);
    const double beta = global_new_r2/global_old_r2;
    handle_boundary_2d(nx, ny, mesh, x, NO_INVERT, PACK);

#if 0
    // Check if the solution has converged
    if(fabs(global_new_r2) < 1.0e-10) {
      global_old_r2 = global_new_r2;
      break;
    }
#endif // if 0

    update_conjugate(nx, ny, beta, r, p);
    handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);

    // Store the old squared residual
    global_old_r2 = global_new_r2;
  }

  *end_niters = ii;
  *end_error = global_old_r2;
#endif // if 0
}

// Initialises the CG solver
double initialise_cg(
    const int nx, const int ny, const double dt, const double conductivity,
    const double heat_capacity, double* p, double* r, const double* x, 
    const double* rho, double* s_x, double* s_y, const double* edgedx, 
    const double* edgedy)
{
  START_PROFILING(&compute_profile);

  // https://inldigitallibrary.inl.gov/sti/3952796.pdf
  // Take the average of the coefficients at the cells surrounding 
  // each face
#pragma omp parallel for
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < (nx+1)-PAD; ++jj) {
      s_x[ind1] = (dt*conductivity*(rho[ind0]+rho[ind0-1]))/
        (2.0*rho[ind0]*rho[ind0-1]*edgedx[jj]*edgedx[jj]*heat_capacity);
    }
  }
#pragma omp parallel for
  for(int ii = PAD; ii < (ny+1)-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      s_y[ind0] = (dt*conductivity*(rho[ind0]+rho[ind0-nx]))/
        (2.0*rho[ind0]*rho[ind0-nx]*edgedy[ii]*edgedy[ii]*heat_capacity);
    }
  }

  double initial_r2 = 0.0;
#pragma omp parallel for reduction(+: initial_r2)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      r[ind0] = x[ind0] -
        ((s_y[ind0]+s_x[ind1]+1.0+s_x[ind1+1]+s_y[ind0+nx])*x[ind0]
         - s_y[ind0]*x[ind0-nx]
         - s_x[ind1]*x[ind0-1] 
         - s_x[ind1+1]*x[ind0+1]
         - s_y[ind0+nx]*x[ind0+nx]);
      p[ind0] = r[ind0];
      initial_r2 += r[ind0]*r[ind0];
    }
  }

  STOP_PROFILING(&compute_profile, "initialise cg");
  return initial_r2;
}

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, const double* s_x, const double* s_y,
    double* p, double* Ap)
{
  START_PROFILING(&compute_profile);

  // You don't need to use a matrix as the band matrix is fully predictable
  // from the 5pt stencil
  double pAp = 0.0;
#pragma omp parallel for reduction(+: pAp)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      Ap[ind0] = 
        (s_y[ind0]+s_x[ind1]+1.0+s_x[ind1+1]+s_y[ind0+nx])*p[ind0]
        - s_y[ind0]*p[ind0-nx]
        - s_x[ind1]*p[ind0-1] 
        - s_x[ind1+1]*p[ind0+1]
        - s_y[ind0+nx]*p[ind0+nx];
      pAp += p[ind0]*Ap[ind0];
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
      x[ind0] += alpha*p[ind0];
      r[ind0] -= alpha*Ap[ind0];
      new_r2 += r[ind0]*r[ind0];
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
      p[ind0] = r[ind0] + beta*p[ind0];
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

