#include "nodes_data.h"
#include "../comms.h"
#include "../mesh.h"
#include "../params.h"
#include "../shared.h"
#include "nodes_interface.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef SILO
#include <silo.h>
#endif

// Build an unstructured mesh
void build_unstructured_quad_mesh(NodesMesh* mesh, const int global_nx,
                                  const int global_ny, const int nx,
                                  const int ny);

/* Implementing a very manual approach to this in order to get something
 * working and then it can be incorporated into the multi-package application
 * at some point */

// Initialises the data specific to the nodes application
void initialise_nodes_data(const int nx, const int ny, NodesData* nodes_data,
                           const char* nodes_params) {
  nodes_data->conductivity = get_double_parameter("conductivity", nodes_params);
  nodes_data->heat_capacity =
      get_double_parameter("heat_capacity", nodes_params);
  allocate_data(&nodes_data->b, nx * ny);
}

// Build a 2d rectilinear mesh of quads
void initialise_rectilinear_quad_mesh_2d(NodesMesh* nmesh,
                                         Mesh* mesh) {
  /* Currently this unstructured mesh generation is specific to quads,
   * going to try to work out what needs to happen to make fully unstructured */

  // Just setting all cells to have same number of edges
  const int nx = mesh->local_nx;
  const int ny = mesh->local_nx;
  const int pad = mesh->pad;
  const int global_nx = mesh->global_nx;
  const int global_ny = mesh->global_nx;
  const double width = mesh->width;
  const double height = mesh->height;

  if (NEDGES != 4) {
    TERMINATE("Only implemented for quads.");
  }

  nmesh->nedges = 2 * nx * ny + nx + ny;

  allocate_data(&nmesh->vertices_x, (nx + 1) * (ny + 1));
  allocate_data(&nmesh->vertices_y, (nx + 1) * (ny + 1));

  // Ordered by edge_vertex_0 is closest to bottom left
  allocate_int_data(&nmesh->edge_vertex0,
                    nmesh->nedges);
  allocate_int_data(&nmesh->edge_vertex1,
                    nmesh->nedges);
  allocate_int_data(&nmesh->edges_cells,
                    nmesh->nedges * NCELLS_PER_EDGE);

  // Construct the list of vertices contiguously, currently Cartesian
  for (int ii = 0; ii < (ny + 1); ++ii) {
    for (int jj = 0; jj < (nx + 1); ++jj) {
      const int index = (ii) * (nx + 1) + (jj);
      nmesh->vertices_x[index] =
          (double)((jj)-pad) * (width / (double)global_nx);
      nmesh->vertices_y[index] =
          (double)((ii)-pad) * (height / (double)global_ny);
    }
  }

  // Calculate the vertices connecting each edge
  for (int ii = 0; ii < ny + 1; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int edge_index = (ii) * (2 * nx + 1) + (jj);
      nmesh->edge_vertex0[edge_index] = (ii) * (nx + 1) + (jj);
      nmesh->edge_vertex1[edge_index] = (ii) * (nx + 1) + (jj + 1);
    }
    if (ii < ny) {
      for (int jj = 0; jj < nx + 1; ++jj) {
        const int edge_index = (ii) * (2 * nx + 1) + (jj) + nx;
        nmesh->edge_vertex0[edge_index] = (ii) * (nx + 1) + (jj);
        nmesh->edge_vertex1[edge_index] =
            (ii + 1) * (nx + 1) + (jj);
      }
    }
  }

  // Calculate the cells connected to edges, as a neighbour list.
  for (int ii = 0; ii < ny + 1; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int edge_index = (ii) * (2 * nx + 1) + (jj);
      nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 0] =
          (ii)*nx + (jj);
      nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 1] =
          (ii - 1) * nx + (jj);
    }
    if (ii < ny) {
      for (int jj = 0; jj < nx + 1; ++jj) {
        const int edge_index = (ii) * (2 * nx + 1) + (jj) + nx;
        nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 0] =
            (ii)*nx + (jj);
        nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 1] =
            (ii)*nx + (jj - 1);
      }
    }
  }

  // TODO: Make sure that the memory order of the cells_edges array is
  // optimal for all architectures
  allocate_int_data(&nmesh->cells_edges, NEDGES * nx * ny);
  allocate_data(&nmesh->cell_centroids_x, nx * ny);
  allocate_data(&nmesh->cell_centroids_y, nx * ny);
  allocate_data(&nmesh->volume, nx * ny);

  // Initialise cells connecting edges
  for (int ii = 0; ii < ny; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      nmesh->cells_edges[(BOTTOM)*nx * ny + (ii)*nx + (jj)] =
          (ii) * (2 * nx + 1) + (jj);
      nmesh->cells_edges[(LEFT)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(BOTTOM)*nx * ny + (ii)*nx + (jj)] +
          nx;
      nmesh->cells_edges[(RIGHT)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(LEFT)*nx * ny + (ii)*nx + (jj)] + 1;
      nmesh->cells_edges[(TOP)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(RIGHT)*nx * ny + (ii)*nx + (jj)] + nx;
    }
  }

  // Find the (x,y) location of each of the cell centroids, and the cell volume
  for (int ii = 0; ii < ny; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii)*nx + (jj);

      double A = 0.0;
      double c_x_factor = 0.0;
      double c_y_factor = 0.0;

      for (int kk = 0; kk < NEDGES; ++kk) {
        int edge_index =
            nmesh->cells_edges[(kk)*nx * ny + (cell_index)];
        int edge_vertex0 = nmesh->edge_vertex0[edge_index];
        int edge_vertex1 = nmesh->edge_vertex1[edge_index];

        // The top and left vertices need to be ordered backwards to ensure
        // correct counter-clockwise access
        double x0 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_x[edge_vertex1]
                        : nmesh->vertices_x[edge_vertex0];
        double y0 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_y[edge_vertex1]
                        : nmesh->vertices_y[edge_vertex0];
        double x1 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_x[edge_vertex0]
                        : nmesh->vertices_x[edge_vertex1];
        double y1 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_y[edge_vertex0]
                        : nmesh->vertices_y[edge_vertex1];

        A += 0.5 * (x0 * y1 - x1 * y0);
        c_x_factor += (x0 + x1) * (x0 * y1 - x1 * y0);
        c_y_factor += (y0 + y1) * (x0 * y1 - x1 * y0);
      }

      // NOTE: This calculation of the volume is actually general to all
      // simple polygons...
      nmesh->volume[(ii)*nx + (jj)] = A;
      nmesh->cell_centroids_x[cell_index] =
          (1.0 / (6.0 * A)) * c_x_factor;
      nmesh->cell_centroids_y[cell_index] =
          (1.0 / (6.0 * A)) * c_y_factor;
    }
  }
}

// Build a 2d curvilinear mesh of quads
void initialise_curvilinear_quad_mesh_2d(NodesMesh* nmesh, Mesh* mesh) {
  /* Currently this unstructured mesh generation is specific to quads,
   * going to try to work out what needs to happen to make fully unstructured */

  // Just setting all cells to have same number of edges
  const int nx = mesh->local_nx;
  const int ny = mesh->local_nx;
  const int pad = mesh->pad;
  const int global_nx = mesh->global_nx;
  const int global_ny = mesh->global_nx;
  const double width = mesh->width;
  const double height = mesh->height;

  if (NEDGES != 4) {
    TERMINATE("Only implemented for quads.");
  }

  nmesh->nedges = 2 * nx * ny + nx + ny;

  allocate_data(&nmesh->vertices_x, (nx + 1) * (ny + 1));
  allocate_data(&nmesh->vertices_y, (nx + 1) * (ny + 1));

  // Ordered by edge_vertex_0 is closest to bottom left
  allocate_int_data(&nmesh->edge_vertex0,
                    nmesh->nedges);
  allocate_int_data(&nmesh->edge_vertex1,
                    nmesh->nedges);
  allocate_int_data(&nmesh->edges_cells,
                    nmesh->nedges * NCELLS_PER_EDGE);

  // Construct the list of vertices contiguously, currently Cartesian
  for (int ii = 0; ii < (ny + 1); ++ii) {
    for (int jj = 0; jj < (nx + 1); ++jj) {
      const int xskip = (jj % 2 == 0);
      const int yskip = (ii % 2 == 0);
      const double cell_width = (width / (double)global_nx);
      const double cell_height = (height / (double)global_ny);
      const int index = (ii) * (nx + 1) + (jj);
      nmesh->vertices_x[index] =
          (double)((jj)-pad) * cell_width +
          (cell_width * (yskip & xskip ? 0.8 : 1.2));
      nmesh->vertices_y[index] =
          (double)((ii)-pad) * cell_height +
          (cell_height * (yskip & xskip ? 0.8 : 1.2));
    }
  }

  // Calculate the vertices connecting each edge
  for (int ii = 0; ii < ny + 1; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int edge_index = (ii) * (2 * nx + 1) + (jj);
      nmesh->edge_vertex0[edge_index] = (ii) * (nx + 1) + (jj);
      nmesh->edge_vertex1[edge_index] = (ii) * (nx + 1) + (jj + 1);
    }
    if (ii < ny) {
      for (int jj = 0; jj < nx + 1; ++jj) {
        const int edge_index = (ii) * (2 * nx + 1) + (jj) + nx;
        nmesh->edge_vertex0[edge_index] = (ii) * (nx + 1) + (jj);
        nmesh->edge_vertex1[edge_index] =
            (ii + 1) * (nx + 1) + (jj);
      }
    }
  }

  // Calculate the cells connected to edges, as a neighbour list.
  for (int ii = 0; ii < ny + 1; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int edge_index = (ii) * (2 * nx + 1) + (jj);
      nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 0] =
          (ii)*nx + (jj);
      nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 1] =
          (ii - 1) * nx + (jj);
    }
    if (ii < ny) {
      for (int jj = 0; jj < nx + 1; ++jj) {
        const int edge_index = (ii) * (2 * nx + 1) + (jj) + nx;
        nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 0] =
            (ii)*nx + (jj);
        nmesh->edges_cells[(edge_index)*NCELLS_PER_EDGE + 1] =
            (ii)*nx + (jj - 1);
      }
    }
  }

  // TODO: Make sure that the memory order of the cells_edges array is
  // optimal for all architectures
  allocate_int_data(&nmesh->cells_edges, NEDGES * nx * ny);
  allocate_data(&nmesh->cell_centroids_x, nx * ny);
  allocate_data(&nmesh->cell_centroids_y, nx * ny);
  allocate_data(&nmesh->volume, nx * ny);

  // Initialise cells connecting edges
  for (int ii = 0; ii < ny; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      nmesh->cells_edges[(BOTTOM)*nx * ny + (ii)*nx + (jj)] =
          (ii) * (2 * nx + 1) + (jj);
      nmesh->cells_edges[(LEFT)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(BOTTOM)*nx * ny + (ii)*nx + (jj)] +
          nx;
      nmesh->cells_edges[(RIGHT)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(LEFT)*nx * ny + (ii)*nx + (jj)] + 1;
      nmesh->cells_edges[(TOP)*nx * ny + (ii)*nx + (jj)] =
          nmesh->cells_edges[(RIGHT)*nx * ny + (ii)*nx + (jj)] + nx;
    }
  }

  // Find the (x,y) location of each of the cell centroids, and the cell volume
  for (int ii = 0; ii < ny; ++ii) {
    for (int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii)*nx + (jj);

      double A = 0.0;
      double c_x_factor = 0.0;
      double c_y_factor = 0.0;

      for (int kk = 0; kk < NEDGES; ++kk) {
        int edge_index =
            nmesh->cells_edges[(kk)*nx * ny + (cell_index)];
        int edge_vertex0 = nmesh->edge_vertex0[edge_index];
        int edge_vertex1 = nmesh->edge_vertex1[edge_index];

        // The top and left vertices need to be ordered backwards to ensure
        // correct counter-clockwise access
        double x0 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_x[edge_vertex1]
                        : nmesh->vertices_x[edge_vertex0];
        double y0 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_y[edge_vertex1]
                        : nmesh->vertices_y[edge_vertex0];
        double x1 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_x[edge_vertex0]
                        : nmesh->vertices_x[edge_vertex1];
        double y1 = (kk == TOP || kk == LEFT)
                        ? nmesh->vertices_y[edge_vertex0]
                        : nmesh->vertices_y[edge_vertex1];

        A += 0.5 * (x0 * y1 - x1 * y0);
        c_x_factor += (x0 + x1) * (x0 * y1 - x1 * y0);
        c_y_factor += (y0 + y1) * (x0 * y1 - x1 * y0);
      }

      // NOTE: This calculation of the volume is actually general to all
      // simple polygons...
      nmesh->volume[(cell_index)] = A;
      nmesh->cell_centroids_x[(cell_index)] =
          (1.0 / (6.0 * A)) * c_x_factor;
      nmesh->cell_centroids_y[(cell_index)] =
          (1.0 / (6.0 * A)) * c_y_factor;
    }
  }
}

// Considering the writing of curvilinear data into a silo file
void write_quad_data_to_visit(const int nx, const int ny, const int step,
                              double* vertices_x, double* vertices_y,
                              const double* data) {
#ifdef ENABLE_SILO
  char filename[MAX_STR_LEN];
  sprintf(filename, "output%04d.silo", step);

  DBfile* dbfile =
      DBCreate(filename, DB_CLOBBER, DB_LOCAL, "simulation time step", DB_HDF5);

  int dims[] = {nx + 1, ny + 1};
  int ndims = 2;
  double* coords[] = {(double*)vertices_x, (double*)vertices_y};
  DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims, DB_DOUBLE,
                DB_NONCOLLINEAR, NULL);

  int dims_nodal[] = {nx + 1, ny + 1};
  DBPutQuadvar1(dbfile, "nodal", "quadmesh", data, dims_nodal, ndims, NULL, 0,
                DB_DOUBLE, DB_ZONECENT, NULL);
  DBClose(dbfile);
#endif
}
