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
    const double* edgedx, const double* edgedy)
{
  // Store initial residual
  double local_old_r2 = initialise_cg(
      nx, ny, dt, heat_capacity, conductivity, p, r, x, rho, s_x, s_y, edgedx, 
      edgedy);
  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, x, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for(ii = 0; ii < max_inners; ++ii) {

    const double local_pAp = calculate_pAp(nx, ny, s_x, s_y, p, Ap);
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
    const double* edgedy)
{
  START_PROFILING(&compute_profile);

  // Going to initialise the coefficients here. This ensures that if the 
  // density or mesh were changed by another package, that the coefficients
  // are updated accordingly, making performance evaluation fairer.
  
  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int cell_index = (ii)*nx+(jj);
      double k = 0.0;
      
      // TODO: In the future nedges will be dependent upon cell.
      for(int ff = 0; ff < nedges; ++ff) {
        const int neighbour_index = neighbours[(ff)*nx*ny+(cell_index)];
        const double density0 = rho[(cell_index)];
        const double density1 = rho[(neighbour_index)];

        // Fetch the vertices for the current edge
        const int edge_index = cells_edges[(ff)*nx*ny+(cell_index)];
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the distance between centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];
        const double cell_dx = (neighbour_centroid_x-cell_centroid_x);
        const double cell_dy = (neighbour_centroid_y-cell_centroid_y);
        const double ds = sqrt(cell_dx*cell_dx+cell_dy*cell_dy);

        // Calculate the unit vector joining cell centroids
        const double es_x = cell_dx/ds;
        const double es_y = cell_dy/ds;

        // Calculate the area vector
        const double A_x = -(vertices_x[vertex1]-vertices_x[vertex0]);
        const double A_y = vertices_y[vertex1]-vertices_y[vertex0];

        const double edge_density = (2.0*density0*density1)/(density0+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

      }
    }
  }

  // https://inldigitallibrary.inl.gov/sti/3952796.pdf
  // Take the average of the coefficients at the cells surrounding 
  // each edge
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
      r[(ii)*nx+(jj)] = x[(ii)*nx+(jj)] -
        ((s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
          s_x[(ii)*(nx+1)+(jj+1)]+s_y[(ii+1)*nx+(jj)])*x[(ii)*nx+(jj)]
         - s_y[(ii)*nx+(jj)]*x[(ii-1)*nx+(jj)]
         - s_x[(ii)*(nx+1)+(jj)]*x[(ii)*nx+(jj-1)] 
         - s_x[(ii)*(nx+1)+(jj+1)]*x[(ii)*nx+(jj+1)]
         - s_y[(ii+1)*nx+(jj)]*x[(ii+1)*nx+(jj)]);
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
      Ap[(ii)*nx+(jj)] = 
        (s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
         s_x[(ii)*(nx+1)+(jj+1)]+s_y[(ii+1)*nx+(jj)])*p[(ii)*nx+(jj)]
        - s_y[(ii)*nx+(jj)]*p[(ii-1)*nx+(jj)]
        - s_x[(ii)*(nx+1)+(jj)]*p[(ii)*nx+(jj-1)] 
        - s_x[(ii)*(nx+1)+(jj+1)]*p[(ii)*nx+(jj+1)]
        - s_y[(ii+1)*nx+(jj)]*p[(ii+1)*nx+(jj)];
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

