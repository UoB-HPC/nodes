#include "nodes.h"
#include "../../../comms.h"
#include "../../../profiler.h"
#include "../../nodes_data.h"
#include "../../nodes_interface.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Solve the unstructured diffusion problem
void solve_unstructured_diffusion_2d(
    const int nx, const int ny, const int pad, Mesh* mesh,
    NodesMesh* nmesh, const int max_inners, const double dt,
    const double heat_capacity, const double conductivity, double* temperature,
    double* b, double* r, double* p, double* rho, double* Ap, int* end_niters,
    double* end_error, double* reduce_array) {

  // Store initial residual
  calculate_rhs(
      nx, ny, pad, heat_capacity, conductivity, dt, nmesh->volume,
      rho, temperature, b, nmesh->edge_vertex0,
      nmesh->edge_vertex1, nmesh->cell_centroids_x,
      nmesh->cell_centroids_y, nmesh->vertices_x,
      nmesh->vertices_y, nmesh->cells_edges,
      nmesh->edges_cells);

  double local_old_r2 = initialise_cg(
      nx, ny, pad, dt, conductivity, heat_capacity, p, r, temperature,
      nmesh->volume, b, rho, nmesh->cells_edges,
      nmesh->edge_vertex0, nmesh->edge_vertex1,
      nmesh->vertices_x, nmesh->vertices_y,
      nmesh->cell_centroids_x, nmesh->cell_centroids_y,
      nmesh->edges_cells);

  double global_old_r2 = reduce_all_sum(local_old_r2);

  handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);
  handle_boundary_2d(nx, ny, mesh, temperature, NO_INVERT, PACK);

  // TODO: Can one of the allreduces be removed with kernel fusion?
  int ii = 0;
  for (ii = 0; ii < max_inners; ++ii) {

    const double local_pAp = calculate_pAp(
        nx, ny, pad, p, Ap, dt, conductivity, heat_capacity, temperature,
        nmesh->volume, rho, nmesh->cells_edges,
        nmesh->edge_vertex0, nmesh->edge_vertex1,
        nmesh->vertices_x, nmesh->vertices_y,
        nmesh->cell_centroids_x,
        nmesh->cell_centroids_y, nmesh->edges_cells);

    const double global_pAp = reduce_all_sum(local_pAp);
    const double alpha = global_old_r2 / global_pAp;

    const double local_new_r2 =
        calculate_new_r2(nx, ny, pad, alpha, temperature, p, r, Ap);
    const double global_new_r2 = reduce_all_sum(local_new_r2);
    const double beta = global_new_r2 / global_old_r2;
    handle_boundary_2d(nx, ny, mesh, temperature, NO_INVERT, PACK);

    // Check if the solution has converged
    if (fabs(global_new_r2) < EPS) {
      global_old_r2 = global_new_r2;
      break;
    }

    update_conjugate(nx, ny, pad, beta, r, p);
    handle_boundary_2d(nx, ny, mesh, p, NO_INVERT, PACK);

    // Store the old squared residual
    global_old_r2 = global_new_r2;
  }

  *end_niters = ii;
  *end_error = global_old_r2;
}

// Calculate the RHS including the unstructured correction term
void calculate_rhs(const int nx, const int ny, const int pad,
                   const double heat_capacity, const double conductivity,
                   const double dt, const double* volume, const double* rho,
                   const double* temperature, double* b,
                   const int* edge_vertex0, const int* edge_vertex1,
                   const double* cell_centroids_x,
                   const double* cell_centroids_y, const double* vertices_x,
                   const double* vertices_y, const int* cells_edges,
                   const int* edges_cells) {
/*
   Note here that the temperature is the guessed temperature.
   */

// Find the RHS that includes the unstructured mesh correction
#pragma omp parallel for
  for (int ii = pad; ii < ny - pad; ++ii) {
#pragma omp simd
    for (int jj = pad; jj < nx - pad; ++jj) {

      // Fetch the cell centered values
      const int cell_index = (ii)*nx + (jj);
      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      // Calculate the cell centroids
      const double cell_centroid_x = cell_centroids_x[(cell_index)];
      const double cell_centroid_y = cell_centroids_y[(cell_index)];

      /*
       *    Performing least squares approximation to get unknown d
       *    d = [ dphi/dx, dphi/dy ]
       *    M = [ (dx0, dx1 ..., dx_ff) (dy0, dy1, ..., dy_ff) ]
       *    del(phi) = [ phi1-phi0, phi2-phi0, ..., phi_ff-phi0 ]
       *    d = (M^T.M)^(-1).(M^T).del(phi)
       */

      // Calculate the coefficents to matrix M
      double MTM[3] = {0.0}; // Describes the three unique quantities in (M^T.M)
      double MT_del_phi[2] = {0.0};
      double coeff[2] = {0.0};

      // Calculate the coefficients for all edges
      for (int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx * ny + (cell_index)];
        const int neighbour_index =
            (edges_cells[edge_index * NCELLS_PER_EDGE] == cell_index)
                ? edges_cells[edge_index * NCELLS_PER_EDGE + 1]
                : edges_cells[edge_index * NCELLS_PER_EDGE];

        // Calculate the vector pointing between the cell centroids
        double es_x = (cell_centroids_x[(neighbour_index)] - cell_centroid_x);
        double es_y = (cell_centroids_y[(neighbour_index)] - cell_centroid_y);
        const double centroid_distance = sqrt(es_x * es_x + es_y * es_y);
        es_x /= centroid_distance;
        es_y /= centroid_distance;

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the area vector, even though vertices aren't ordered well
        double A_x = (vertices_y[vertex1] - vertices_y[vertex0]);
        double A_y = -(vertices_x[vertex1] - vertices_x[vertex0]);
        if ((A_x * es_x + A_y * es_y) < 0.0) {
          A_x = -A_x;
          A_y = -A_y;
        }

        // Calculate the gradient matrix
        const double phi0 = temperature[(cell_index)];
        const double phi_ff = temperature[(neighbour_index)];
        MTM[0] += es_x * es_x;
        MTM[1] += es_x * es_y;
        MTM[2] += es_y * es_y;
        MT_del_phi[0] += es_x * (phi_ff - phi0);
        MT_del_phi[1] += es_y * (phi_ff - phi0);

        // Calculate the coefficients of transformed shape
        const double density1 = rho[(neighbour_index)];
        const double edge_density =
            (2.0 * density * density1) / (density + density1);
        const double diffusion_coeff =
            conductivity / (edge_density * heat_capacity);
        const double gam = (A_x * A_x + A_y * A_y) / (A_x * es_x + A_y * es_y);
        coeff[0] += diffusion_coeff * (A_x - es_x * gam);
        coeff[1] += diffusion_coeff * (A_y - es_y * gam);
      }

      // Solve the equation for the temperature gradients
      const double MTM_det = (1.0 / (MTM[0] * MTM[2] - MTM[1] * MTM[1]));
      const double temp_grad_cell_x =
          MTM_det * (MT_del_phi[0] * MTM[2] - MT_del_phi[1] * MTM[1]);
      const double temp_grad_cell_y =
          MTM_det * (MT_del_phi[1] * MTM[0] - MT_del_phi[0] * MTM[1]);

      // TODO: SHOULD THERE BE A COEFFICIENT FOR TAU?
      const double tau =
          temp_grad_cell_x * coeff[0] + temp_grad_cell_y * coeff[1];
      b[(cell_index)] = temperature[(cell_index)] + (dt / (density * V)) * tau;
    }
  }
}

// Initialises the CG solver
double initialise_cg(const int nx, const int ny, const int pad, const double dt,
                     const double conductivity, const double heat_capacity,
                     double* p, double* r, const double* temperature,
                     const double* volume, const double* b, const double* rho,
                     const int* cells_edges, const int* edge_vertex0,
                     const int* edge_vertex1, const double* vertices_x,
                     const double* vertices_y, const double* cell_centroids_x,
                     const double* cell_centroids_y, const int* edges_cells) {
  START_PROFILING(&compute_profile);

  // Going to initialise the coefficients here. This ensures that if the
  // density or mesh were changed by another package, that the coefficients
  // are updated accordingly, making performance evaluation fairer.

  double initial_r2 = 0.0;
#pragma omp parallel for reduction(+ : initial_r2)
  for (int ii = pad; ii < ny - pad; ++ii) {
#pragma omp simd
    for (int jj = pad; jj < nx - pad; ++jj) {
      const int cell_index = (ii)*nx + (jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      // Calculate the cell centroids
      const double cell_centroid_x = cell_centroids_x[(cell_index)];
      const double cell_centroid_y = cell_centroids_y[(cell_index)];

      double neighbour_coeff_total = 0.0;
      double neighbour_contribution = 0.0;

      for (int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx * ny + (cell_index)];
        const int neighbour_index =
            (edges_cells[edge_index * NCELLS_PER_EDGE] == cell_index)
                ? edges_cells[edge_index * NCELLS_PER_EDGE + 1]
                : edges_cells[edge_index * NCELLS_PER_EDGE];

        // Calculate the unit vector pointing between the cell centroids
        double es_x = (cell_centroids_x[(neighbour_index)] - cell_centroid_x);
        double es_y = (cell_centroids_y[(neighbour_index)] - cell_centroid_y);
        const double centroid_distance = sqrt(es_x * es_x + es_y * es_y);
        es_x /= centroid_distance;
        es_y /= centroid_distance;

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the area vector, even though vertices aren't ordered well
        double A_x = (vertices_y[vertex1] - vertices_y[vertex0]);
        double A_y = -(vertices_x[vertex1] - vertices_x[vertex0]);
        if ((A_x * es_x + A_y * es_y) < 0.0) {
          A_x = -A_x;
          A_y = -A_y;
        }

        // Calculate the diffusion coefficient
        const double edge_density = (2.0 * density * rho[(neighbour_index)]) /
                                    (density + rho[(neighbour_index)]);
        const double diffusion_coeff =
            conductivity / (edge_density * heat_capacity);
        const double neighbour_coeff =
            (dt * diffusion_coeff * (A_x * A_x + A_y * A_y)) /
            (V * centroid_distance * (A_x * es_x + A_y * es_y));
        neighbour_contribution +=
            temperature[(neighbour_index)] * neighbour_coeff;
        neighbour_coeff_total += neighbour_coeff;
      }

      r[(cell_index)] = b[(cell_index)] - ((neighbour_coeff_total + 1.0) *
                                               temperature[(cell_index)] -
                                           neighbour_contribution);
      p[(cell_index)] = r[(cell_index)];
      initial_r2 += r[(cell_index)] * r[(cell_index)];
    }
  }

  STOP_PROFILING(&compute_profile, "initialise cg");
  return initial_r2;
}

// Calculates a value for alpha
double calculate_pAp(const int nx, const int ny, const int pad, double* p,
                     double* Ap, const double dt, const double conductivity,
                     const double heat_capacity, const double* temperature,
                     const double* volume, const double* rho,
                     const int* cells_edges, const int* edge_vertex0,
                     const int* edge_vertex1, const double* vertices_x,
                     const double* vertices_y, const double* cell_centroids_x,
                     const double* cell_centroids_y, const int* edges_cells) {
  START_PROFILING(&compute_profile);

  double pAp = 0.0;
#pragma omp parallel for reduction(+ : pAp)
  for (int ii = pad; ii < ny - pad; ++ii) {
#pragma omp simd
    for (int jj = pad; jj < nx - pad; ++jj) {
      const int cell_index = (ii)*nx + (jj);

      const double density = rho[(cell_index)];
      const double V = volume[(cell_index)];

      // Calculate the cell centroids
      const double cell_centroid_x = cell_centroids_x[(cell_index)];
      const double cell_centroid_y = cell_centroids_y[(cell_index)];

      double neighbour_coeff_total = 0.0;
      double neighbour_contribution = 0.0;

      for (int ee = 0; ee < NEDGES; ++ee) {
        const int edge_index = cells_edges[(ee)*nx * ny + (cell_index)];
        const int neighbour_index =
            (edges_cells[edge_index * NCELLS_PER_EDGE] == cell_index)
                ? edges_cells[edge_index * NCELLS_PER_EDGE + 1]
                : edges_cells[edge_index * NCELLS_PER_EDGE];

        // Calculate the unit vector pointing between the cell centroids
        double es_x = (cell_centroids_x[(neighbour_index)] - cell_centroid_x);
        double es_y = (cell_centroids_y[(neighbour_index)] - cell_centroid_y);
        const double centroid_distance = sqrt(es_x * es_x + es_y * es_y);
        es_x /= centroid_distance;
        es_y /= centroid_distance;

        // Calculate the edge differentials
        const int vertex0 = edge_vertex0[(edge_index)];
        const int vertex1 = edge_vertex1[(edge_index)];

        // Calculate the area vector, even though vertices aren't ordered well
        double A_x = (vertices_y[vertex1] - vertices_y[vertex0]);
        double A_y = -(vertices_x[vertex1] - vertices_x[vertex0]);
        if ((A_x * es_x + A_y * es_y) < 0.0) {
          A_x = -A_x;
          A_y = -A_y;
        }

        // Calculate the diffusion coefficient
        const double edge_density = (2.0 * density * rho[(neighbour_index)]) /
                                    (density + rho[(neighbour_index)]);
        const double diffusion_coeff =
            conductivity / (edge_density * heat_capacity);
        const double neighbour_coeff =
            (dt * diffusion_coeff * (A_x * A_x + A_y * A_y)) /
            (V * centroid_distance * (A_x * es_x + A_y * es_y));
        neighbour_contribution += p[(neighbour_index)] * neighbour_coeff;
        neighbour_coeff_total += neighbour_coeff;
      }

      Ap[(cell_index)] = ((neighbour_coeff_total + 1.0) * p[(cell_index)] -
                          neighbour_contribution);
      pAp += p[(cell_index)] * Ap[(cell_index)];
    }
  }

  STOP_PROFILING(&compute_profile, "calculate alpha");
  return pAp;
}

// Updates the current guess using the calculated alpha
double calculate_new_r2(const int nx, const int ny, const int pad, double alpha,
                        double* temperature, double* p, double* r, double* Ap) {
  START_PROFILING(&compute_profile);

  double new_r2 = 0.0;

  for (int ii = pad; ii < ny - pad; ++ii) {
    for (int jj = pad; jj < nx - pad; ++jj) {
      temperature[(ii)*nx + (jj)] += alpha * p[(ii)*nx + (jj)];
      r[(ii)*nx + (jj)] -= alpha * Ap[(ii)*nx + (jj)];
      new_r2 += r[(ii)*nx + (jj)] * r[(ii)*nx + (jj)];
    }
  }

  STOP_PROFILING(&compute_profile, "calculate new r2");
  return new_r2;
}

// Updates the conjugate from the calculated beta and residual
void update_conjugate(const int nx, const int ny, const int pad,
                      const double beta, const double* r, double* p) {
  START_PROFILING(&compute_profile);
  for (int ii = pad; ii < ny - pad; ++ii) {
    for (int jj = pad; jj < nx - pad; ++jj) {
      p[(ii)*nx + (jj)] = r[(ii)*nx + (jj)] + beta * p[(ii)*nx + (jj)];
    }
  }
  STOP_PROFILING(&compute_profile, "update conjugate");
}
