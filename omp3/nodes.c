#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nodes.h"
#include "../nodes_data.h"
#include "../nodes_interface.h"
#include "../../profiler.h"
#include "../../comms.h"

// Solve the unstructured diffusion problem
void solve_unstructured_diffusion_2d(
    const int nx, const int ny, Mesh* mesh, UnstructuredMesh* unstructured_mesh, 
    const int max_inners, const double dt, const double heat_capacity, 
    const double conductivity, double* temperature, double* b, double* r, double* p, 
    double* rho, double* Ap, int* end_niters, double* end_error, double* reduce_array)
{
  // Store initial residual
  calculate_rhs( 
      nx, ny, heat_capacity, conductivity, dt, 
      unstructured_mesh->volume, rho, temperature, b, 
      unstructured_mesh->edge_vertex0, unstructured_mesh->edge_vertex1, 
      unstructured_mesh->cell_centroids_x, unstructured_mesh->cell_centroids_y, 
      unstructured_mesh->vertices_x, unstructured_mesh->vertices_y, 
      unstructured_mesh->cells_edges, unstructured_mesh->edges_cells);

  write_quad_data_to_visit(
      mesh->local_nx, mesh->local_ny, 1, unstructured_mesh->vertices_x, 
      unstructured_mesh->vertices_y, b);

  double local_old_r2 = initialise_cg(
      nx, ny, dt, conductivity, heat_capacity, 
      p, r, temperature, unstructured_mesh->volume, b, 
      rho, unstructured_mesh->cells_edges, unstructured_mesh->edge_vertex0, 
      unstructured_mesh->edge_vertex1, unstructured_mesh->vertices_x, 
      unstructured_mesh->vertices_y, unstructured_mesh->cell_centroids_x, 
      unstructured_mesh->cell_centroids_y, unstructured_mesh->edges_cells);

  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, temperature, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for(ii = 0; ii < max_inners; ++ii) {

    const double local_pAp = calculate_pAp(
        nx, ny, p, Ap, dt, conductivity, heat_capacity, temperature, 
        unstructured_mesh->volume, rho, unstructured_mesh->cells_edges, 
        unstructured_mesh->edge_vertex0, unstructured_mesh->edge_vertex1, 
        unstructured_mesh->vertices_x, unstructured_mesh->vertices_y, 
        unstructured_mesh->cell_centroids_x, unstructured_mesh->cell_centroids_y, 
        unstructured_mesh->edges_cells);

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

// Calculate the RHS including the unstructured correction term
void calculate_rhs(
    const int nx, const int ny, const double heat_capacity, const double conductivity, 
    const double dt, const double* volume, const double* rho, const double* temperature, 
    double* b, const int* edge_vertex0, const int* edge_vertex1, 
    const double* cell_centroids_x, const double* cell_centroids_y, const double* vertices_x, 
    const double* vertices_y, const int* cells_edges, const int* edges_cells)
{
  /*
     Note here that the temperature is the guessed temperature.
     */

  // Find the RHS that includes the unstructured mesh correction
  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {

      // Fetch the cell centered values
      const int cell_index = (ii)*nx+(jj);
      double V = volume[(cell_index)];
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
      for(int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx*ny+(cell_index)];
        const int neighbour_index = 
          (edges_cells[edge_index*NCELLS_PER_EDGE] == cell_index) ?
          edges_cells[edge_index*NCELLS_PER_EDGE+1] :
          edges_cells[edge_index*NCELLS_PER_EDGE];

        // Calculate the cell centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the vector pointing between the cell centroids
        const double es_x = (neighbour_centroid_x-cell_centroid_x);
        const double es_y = (neighbour_centroid_y-cell_centroid_y);

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the area vector
        double A_x = (vertices_y[vertex1]-vertices_y[vertex0]);
        double A_y = -(vertices_x[vertex1]-vertices_x[vertex0]);

        if(samesign(es_x , A_y)) {
          A_x = -A_x;
        }
        else if(!samesign(es_y , A_x)) {
          A_y = -A_y;
        }

        // Calculate the gradient matrix
        const double phi0 = temperature[(cell_index)];
        const double phi_ff = temperature[(neighbour_index)];
        MTM[0] += es_x*es_x;
        MTM[1] += es_x*es_y;
        MTM[2] += es_y*es_y;
        MT_del_phi[0] += es_x*(phi_ff-phi0);
        MT_del_phi[1] += es_y*(phi_ff-phi0);

        // Calculate the coefficients of transformed shape
        const double density1 = rho[(neighbour_index)];
        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);
        const double gam = (A_x*A_x+A_y*A_y)/(A_x*es_x+A_y*es_y);
        coeff[0] += diffusion_coeff*(A_x-es_x*gam);
        coeff[1] += diffusion_coeff*(A_y-es_y*gam);
      }

      // Solve the equation for the temperature gradients
      const double MTM_det = (1.0/(MTM[0]*MTM[2]-MTM[1]*MTM[1]));
      const double temp_grad_cell_x = 
        MTM_det*(MT_del_phi[0]*MTM[2]-MT_del_phi[1]*MTM[1]);
      const double temp_grad_cell_y = 
        MTM_det*(MT_del_phi[1]*MTM[0]-MT_del_phi[0]*MTM[1]);

      const double tau = temp_grad_cell_x*coeff[0]+temp_grad_cell_y*coeff[1];
      b[(cell_index)] = (V*density/dt)*temperature[(cell_index)] + tau;
    }
  }
}


// Initialises the CG solver
double initialise_cg(
    const int nx, const int ny, const double dt, const double conductivity,
    const double heat_capacity, double* p, double* r, const double* temperature, 
    const double* volume, const double* b, const double* rho, const int* cells_edges, 
    const int* edge_vertex0, const int* edge_vertex1, const double* vertices_x, 
    const double* vertices_y, const double* cell_centroids_x, 
    const double* cell_centroids_y, const int* edges_cells)
{
  START_PROFILING(&compute_profile);

  // Going to initialise the coefficients here. This ensures that if the 
  // density or mesh were changed by another package, that the coefficients
  // are updated accordingly, making performance evaluation fairer.

  double initial_r2 = 0.0;
  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int cell_index = (ii)*nx+(jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      double neighbour_coeff_total = 0.0;
      double neighbour_contribution = 0.0;

      for(int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx*ny+(cell_index)];
        const int neighbour_index = 
          (edges_cells[edge_index*NCELLS_PER_EDGE] == cell_index) ?
          edges_cells[edge_index*NCELLS_PER_EDGE+1] :
          edges_cells[edge_index*NCELLS_PER_EDGE];

        // Calculate the cell centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the vector pointing between the cell centroids
        const double es_x = (neighbour_centroid_x-cell_centroid_x);
        const double es_y = (neighbour_centroid_y-cell_centroid_y);

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the distance between the centroids
        const double centroid_distance = sqrt(es_x*es_x+es_y*es_y);

        // Calculate the area vector
        const double A_x = (vertices_y[vertex1]-vertices_y[vertex0]);
        const double A_y = -(vertices_x[vertex1]-vertices_x[vertex0]);

        // Calculate the diffusion coefficient
        const double density1 = rho[(neighbour_index)];
        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

        const double neighbour_coeff = 
          (diffusion_coeff*(A_x*A_x+A_y*A_y))/
          (centroid_distance*fabs(A_x*es_x+A_y*es_y));
        neighbour_contribution += temperature[(neighbour_index)]*neighbour_coeff;
        neighbour_coeff_total += neighbour_coeff;
      }

      r[(cell_index)] = b[(cell_index)] - 
        (neighbour_coeff_total+(density*V/dt))*temperature[(cell_index)] - 
        neighbour_contribution;
      p[(cell_index)] = r[(cell_index)];
      initial_r2 += r[(cell_index)]*r[(cell_index)];
    }
  }

  STOP_PROFILING(&compute_profile, "initialise cg");
  return initial_r2;
}

// Calculates a value for alpha
double calculate_pAp(
    const int nx, const int ny, double* p, double* Ap, const double dt, 
    const double conductivity, const double heat_capacity, 
    const double* temperature, const double* volume, const double* rho, 
    const int* cells_edges, const int* edge_vertex0, const int* edge_vertex1, 
    const double* vertices_x, const double* vertices_y,
    const double* cell_centroids_x, const double* cell_centroids_y, 
    const int* edges_cells)
{
  START_PROFILING(&compute_profile);

  double pAp = 0.0;
  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      const int cell_index = (ii)*nx+(jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      double neighbour_coeff_total = 0.0;
      double neighbour_contribution = 0.0;

      for(int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx*ny+(cell_index)];
        const int neighbour_index = 
          (edges_cells[edge_index*NCELLS_PER_EDGE] == cell_index) ?
          edges_cells[edge_index*NCELLS_PER_EDGE+1] :
          edges_cells[edge_index*NCELLS_PER_EDGE];

        // Calculate the cell centroids
        const double cell_centroid_x = cell_centroids_x[(cell_index)];
        const double cell_centroid_y = cell_centroids_y[(cell_index)];
        const double neighbour_centroid_x = cell_centroids_x[(neighbour_index)];
        const double neighbour_centroid_y = cell_centroids_y[(neighbour_index)];

        // Calculate the vector pointing between the cell centroids
        const double es_x = (neighbour_centroid_x-cell_centroid_x);
        const double es_y = (neighbour_centroid_y-cell_centroid_y);

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the distance between the centroids
        const double centroid_distance = sqrt(es_x*es_x+es_y*es_y);

        // Calculate the area vector
        const double A_x = (vertices_y[vertex1]-vertices_y[vertex0]);
        const double A_y = -(vertices_x[vertex1]-vertices_x[vertex0]);

        // Calculate the diffusion coefficient
        const double density1 = rho[(neighbour_index)];
        const double edge_density = (2.0*density*density1)/(density+density1);
        const double diffusion_coeff = conductivity/(edge_density*heat_capacity);

        const double neighbour_coeff = 
          (diffusion_coeff*(A_x*A_x+A_y*A_y))/
          (centroid_distance*fabs(A_x*es_x+A_y*es_y));
        neighbour_contribution += p[(neighbour_index)]*neighbour_coeff;
        neighbour_coeff_total += neighbour_coeff;
      }

      Ap[(cell_index)] = 
        (neighbour_coeff_total+(density*V/dt))*p[(cell_index)]-neighbour_contribution;
      pAp += p[(cell_index)]*Ap[(cell_index)];
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

  for(int ii = PAD; ii < ny-PAD; ++ii) {
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
  for(int ii = PAD; ii < ny-PAD; ++ii) {
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

