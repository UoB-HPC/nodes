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
    double* temperature, double* r, double* p, double* rho, double* s_x, double* s_y, 
    double* Ap, int* end_niters, double* end_error, double* reduce_array,
    const double* edgedx, const double* edgedy)
{
  // Store initial residual
  double local_old_r2 = initialise_cg(
      nx, ny, dt, heat_capacity, conductivity, p, r, temperature, rho, s_x, s_y, edgedx, 
      edgedy);
  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, temperature, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for(ii = 0; ii < max_inners; ++ii) {

    const double local_pAp = calculate_pAp(nx, ny, s_x, s_y, p, Ap);
    const double global_pAp = reduce_all_sum(local_pAp);
    const double alpha = global_old_r2/global_pAp;

    const double local_new_r2 = calculate_new_r2(nx, ny, alpha, temperature, p, r, Ap);
    const double global_new_r2 = reduce_all_sum(local_new_r2);
    const double beta = global_new_r2/global_old_r2;
    handle_boundary_2d(nx, ny, mesh, temperature, NO_INVERT, PACK);

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
    const double heat_capacity, const int nedges, double* p, double* r, 
    const double* temperature, const double* volume, const double* b,
    const double* rho, double* s_x, double* s_y, const double* edgedx, 
    const double* edgedy, const int* cells_edges, const int* edge_vertex0,
    const int* edge_vertex1, const int* vertices_x, const int* vertices_y,
    const int* cell_centroids_x, const int* cell_centroids_y, 
    const int* neighbours)
{
  START_PROFILING(&compute_profile);

  // Going to initialise the coefficients here. This ensures that if the 
  // density or mesh were changed by another package, that the coefficients
  // are updated accordingly, making performance evaluation fairer.
  
  double initial_r2 = 0.0;
#pragma omp parallel for reduction(+: initial_r2)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int cell_index = (ii)*nx+(jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      double cell_coeff = 1;
      double neighbour_contribution = 0;

      for(int ff = 0; ff < nedges; ++ff) {
        const int neighbour_index = neighbours[(ff)*nx*ny+(cell_index)];

        // Calculate the distance between centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the cell differentials
        const double cell_dx = (neighbour_centroid_x-cell_centroid_x);
        const double cell_dy = (neighbour_centroid_y-cell_centroid_y);

        // Calculate the edge differentials
        const int edge_index = cells_edges[(ff)*nx*ny+(cell_index)];
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];
        const double edge_dx = (vertices_x[vertex1]-vertices_x[vertex0]);
        const double edge_dy = (vertices_y[vertex1]-vertices_y[vertex0]);

        // Calculate the centroid distance and length of edge
        const double centroid_distance = sqrt(cell_dx*cell_dx+cell_dy*cell_dy);

        // Calculate the unit vector joining cell centroids
        const double es_x = cell_dx/centroid_distance;
        const double es_y = cell_dy/centroid_distance;

        // Calculate the area vector
        const double A_x = edge_dy;
        const double A_y = -edge_dx;

        // Calculate the diffusion coefficient
        const double density1 = rho[(neighbour_index)];
        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

        const double local_coeff = (dt*diffusion_coeff*(A_x*A_x+A_y*A_y))/
          (centroid_distance*(A_x*es_x+A_y*es_y)*density*V);
        neighbour_contribution = temperature[(neighbour_index)]*local_coeff;
        cell_coeff -= local_coeff;
      }

      r[(ii)*nx+(jj)] = b[(ii)*nx+(jj)] - 
        (neighbour_contribution + cell_coeff*temperature[(cell_index)]);
      p[(ii)*nx+(jj)] = r[(ii)*nx+(jj)];
      initial_r2 += r[(ii)*nx+(jj)]*r[(ii)*nx+(jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "initialise cg");
  return initial_r2;
}

// Calculate the RHS including the unstructured correction term
void calculate_rhs(
    const int nx, const int ny, const int nedges, const double heat_capacity,
    const double conductivity, const double dt, const double* volume, 
    const double* rho, const double* temperature, double* b, 
    const int* edge_vertex0, const int* edge_vertex1,
    const int* cell_centroids_x, const int* cell_centroids_y,
    const int* vertices_x, const int* vertices_y,
    const int* cells_edges, const int* neighbours)
{
  /*
     Note here that the temperature is the guessed temperature.
   */

  // Find the RHS that includes the unstructured mesh correction
  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {

      // Fetch the cell centered values
      const int cell_index = (ii)*nx+(jj);
      double V = volume[cell_index];
      const double density = rho[(cell_index)];

      /*
       *    Performing least squares approximation to get unknown d
       *    d = [ dphi/dx, dphi/dy ]
       *    M = [ (dx0, dx1 ..., dx_ff) (dy0, dy1, ..., dy_ff) ]
       *    del(phi) = [ phi1-phi0, phi2-phi0, ..., phi_ff-phi0 ]
       *    d = (M^T.M)^(-1).(M^T).del(phi)
       */

      // Calculate the coefficents to matrix M
      double MTM[3] = { 0.0 }; // Describes the three unique quantities in (M^T.M)
      double MT_del_phi[2] = { 0.0 };
      double coeff[2] = { 0.0 };

      // Calculate the coefficients for all edges
      for(int ff = 0; ff < nedges; ++ff) {
        const int neighbour_index = neighbours[(ff)*nx*ny+(cell_index)];

        // Calculate the distance between centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the cell differentials
        const double cell_dx = (neighbour_centroid_x-cell_centroid_x);
        const double cell_dy = (neighbour_centroid_y-cell_centroid_y);

        const double phi0 = temperature[(cell_index)];
        const double phi_ff = temperature[(neighbour_index)];
        MTM[0] += cell_dx*cell_dx;
        MTM[1] += cell_dx*cell_dy;
        MTM[2] += cell_dy*cell_dy;
        MT_del_phi[0] += cell_dx*(phi_ff-phi0);
        MT_del_phi[1] += cell_dy*(phi_ff-phi0);

        // Calculate the edge differentials
        const int edge_index = cells_edges[(ff)*nx*ny+(cell_index)];
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];
        const double edge_dx = (vertices_x[vertex1]-vertices_x[vertex0]);
        const double edge_dy = (vertices_y[vertex1]-vertices_y[vertex0]);

        // Calculate the centroid distance and length of edge
        const double centroid_distance = sqrt(cell_dx*cell_dx+cell_dy*cell_dy);

        // Calculate the unit vector joining cell centroids
        const double es_x = cell_dx/centroid_distance;
        const double es_y = cell_dy/centroid_distance;

        // Calculate the area vector
        const double A_x = edge_dy;
        const double A_y = -edge_dx;

        const double density1 = rho[(neighbour_index)];

        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

        coeff[0] = diffusion_coeff*(A_x-(A_x*A_x)/(A_x*es_x));
        coeff[1] = diffusion_coeff*(A_y-(A_y*A_y)/(A_y*es_y));
      }

      // Solve the equation for the temperature gradients
      const double MTM_det = (1.0/(MTM[0]*MTM[2]-MTM[1]*MTM[1]));
      const double temp_grad_cell_x = 
        MTM_det*(MT_del_phi[0]*MTM[2]-MT_del_phi[1]*MTM[1]);
      const double temp_grad_cell_y = 
        MTM_det*(MT_del_phi[1]*MTM[0]-MT_del_phi[0]*MTM[1]);

      const double tau = temp_grad_cell_x*coeff[0]+temp_grad_cell_y*coeff[1];
      b[(cell_index)] = temperature[(cell_index)] - (dt*tau)/(V*density);
    }
  }
}

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, double* p, double* Ap, const double dt, 
    const double conductivity, const double heat_capacity, const int nedges, 
    double* r, const double* temperature, const double* volume,
    const double* rho, double* s_x, double* s_y, const double* edgedx, 
    const double* edgedy, const int* cells_edges, const int* edge_vertex0,
    const int* edge_vertex1, const int* vertices_x, const int* vertices_y,
    const int* cell_centroids_x, const int* cell_centroids_y, 
    const int* neighbours)
{
  START_PROFILING(&compute_profile);

  double pAp = 0.0;
#pragma omp parallel for reduction(+: pAp)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int cell_index = (ii)*nx+(jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      double cell_coeff = 1;
      double neighbour_contribution = 0;

      for(int ff = 0; ff < nedges; ++ff) {
        const int neighbour_index = neighbours[(ff)*nx*ny+(cell_index)];

        // Calculate the distance between centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the cell differentials
        const double cell_dx = (neighbour_centroid_x-cell_centroid_x);
        const double cell_dy = (neighbour_centroid_y-cell_centroid_y);

        // Calculate the edge differentials
        const int edge_index = cells_edges[(ff)*nx*ny+(cell_index)];
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];
        const double edge_dx = (vertices_x[vertex1]-vertices_x[vertex0]);
        const double edge_dy = (vertices_y[vertex1]-vertices_y[vertex0]);

        // Calculate the centroid distance and length of edge
        const double centroid_distance = sqrt(cell_dx*cell_dx+cell_dy*cell_dy);

        // Calculate the unit vector joining cell centroids
        const double es_x = cell_dx/centroid_distance;
        const double es_y = cell_dy/centroid_distance;

        // Calculate the area vector
        const double A_x = edge_dy;
        const double A_y = -edge_dx;

        // Calculate the diffusion coefficient
        const double density1 = rho[(neighbour_index)];
        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

        const double local_coeff = (dt*diffusion_coeff*(A_x*A_x+A_y*A_y))/
          (centroid_distance*(A_x*es_x+A_y*es_y)*density*V);
        neighbour_contribution = temperature[(neighbour_index)]*local_coeff;
        cell_coeff -= local_coeff;
      }

      Ap[(ii)*nx+(jj)] = 
        (neighbour_contribution+cell_coeff*temperature[(cell_index)]);
      pAp += p[(ii)*nx+(jj)]*Ap[(ii)*nx+(jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "calculate alpha");
  return pAp;
}

// Updates the current guess using the calculated alpha
double calculate_new_r2(
    int nx, int ny, double alpha, double* temperature, double* p, double* r, double* Ap)
{
  START_PROFILING(&compute_profile);

  double new_r2 = 0.0;

#pragma omp parallel for reduction(+: new_r2)
  for(int ii = PAD; ii < ny-PAD; ++ii) {
#pragma omp simd
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      temperature[(ii)*nx+(jj)] += alpha*p[(ii)*nx+(jj)];
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

