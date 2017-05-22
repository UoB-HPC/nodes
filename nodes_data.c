#include <stdio.h>
#include <stdlib.h>
#include "nodes_data.h"
#include "nodes_interface.h"
#include "../params.h"
#include "../mesh.h"
#include "../shared.h"
#include "../comms.h"

// Build an unstructured mesh
void build_unstructured_quad_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int nx, const int ny);

/* Implementing a very manual approach to this in order to get something 
 * working and then it can be incorporated into the multi-package application
 * at some point */

void initialise_nodes_data(
    NodesData* nodes_data, const int global_nx, const int global_ny, 
    const int nx, const int ny, const int nneighbours, const char* nodes_params)
{
  nodes_data->conductivity = get_double_parameter("conductivity", nodes_params);
  nodes_data->heat_capacity = get_double_parameter("heat_capacity", nodes_params);
  nodes_data->nneighbours = nneighbours;

  // Describe the neighbours for all real cells
  allocate_int_data(&nodes_data->neighbours_ii, nneighbours*nx*ny);
  allocate_int_data(&nodes_data->neighbours_jj, nneighbours*nx*ny);

  UnstructuredMesh mesh;
  mesh.width = get_double_parameter("width", ARCH_ROOT_PARAMS);
  mesh.height = get_double_parameter("height", ARCH_ROOT_PARAMS);
  build_unstructured_quad_mesh(
      &mesh, global_nx, global_ny, nx, ny);
}

// Build an unstructured mesh
void build_unstructured_quad_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int nx, const int ny)
{
  // TODO: Generalise this code for all types of grid? Or should it be specified
  // to some extent?
  if(NEDGES != 4) {
    TERMINATE("Only implemented for quads currently.");
  }

  mesh->nedges = 2*nx*ny+nx+ny;

  allocate_data(&mesh->vertices_x, (nx+1)*(ny+1));
  allocate_data(&mesh->vertices_y, (nx+1)*(ny+1));

  // Ordered by edge_vertex_0 is closest to bottom left
  allocate_int_data(&mesh->edge_vertex0, mesh->nedges);
  allocate_int_data(&mesh->edge_vertex1, mesh->nedges);

  // Just setting all cells to have same number of edges
  const int nedges = NEDGES;

  // Construct the list of vertices contiguously, currently Cartesian
  for(int ii = 0; ii < (ny+1); ++ii) {
    for(int jj = 0; jj < (nx+1); ++jj) {
      const int index = (ii)*(nx+1)+(jj);
      mesh->vertices_x[index] = (double)((jj)-PAD)*(mesh->width/(double)global_nx);
      mesh->vertices_y[index] = (double)((ii)-PAD)*(mesh->height/(double)global_ny);

      printf("vertex %d (%.4f %.4f)\n", 
          index, mesh->vertices_x[index], mesh->vertices_y[index]);
    }
  }

  // Calculate the vertices connecting each edge, we step through from bottom 
  // left all the way to top right
  for(int ii = 0; ii < ny+1; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int edge_index = (ii)*(2*nx+1)+(jj);
      mesh->edge_vertex0[edge_index] = (ii)*(nx+1)+(jj);
      mesh->edge_vertex1[edge_index] = (ii)*(nx+1)+(jj+1);
    }
    if(ii < ny) {
      for(int jj = 0; jj < nx+1; ++jj) {
        const int edge_index = (ii)*(2*nx+1)+(jj)+nx;
        mesh->edge_vertex0[edge_index] = (ii)*(nx+1)+(jj);
        mesh->edge_vertex1[edge_index] = (ii+1)*(nx+1)+(jj);
      }
    }
  }

  // TODO: Make sure that the memory order of the cells_edges array is
  // optimal for all architectures
  allocate_int_data(&mesh->cells_edges, nedges*nx*ny);

  // Initialise cells connecting edges
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      mesh->cells_edges[BOTTOM*nx*ny+(ii)*nx+(jj)] = (ii)*(2*nx+1)+(jj);
      mesh->cells_edges[LEFT*nx*ny+(ii)*nx+(jj)] = mesh->cells_edges[BOTTOM*nx*ny+(ii)*nx+(jj)]+nx;
      mesh->cells_edges[RIGHT*nx*ny+(ii)*nx+(jj)] = mesh->cells_edges[LEFT*nx*ny+(ii)*nx+(jj)]+1;
      mesh->cells_edges[TOP*nx*ny+(ii)*nx+(jj)] = mesh->cells_edges[RIGHT*nx*ny+(ii)*nx+(jj)]+nx;
    }
  }

  allocate_data(&mesh->cell_centers_x, nx*ny);
  allocate_data(&mesh->cell_centers_y, nx*ny);

  // Find the (x,y) location of each of the cell centers
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii)*nx+(jj);

      double A = 0.0;
      double c_x_factor = 0.0;
      double c_y_factor = 0.0;

      for(int kk = 0; kk < nedges; ++kk) {
        const int edge_index = mesh->cells_edges[(kk)*nx*ny+cell_index];
        const int edge_vertex0 = mesh->edge_vertex0[edge_index];
        const int edge_vertex1 = mesh->edge_vertex1[edge_index];
        const double x0 = mesh->vertices_x[edge_vertex0];
        const double y0 = mesh->vertices_y[edge_vertex0];
        const double x1 = mesh->vertices_x[edge_vertex1];
        const double y1 = mesh->vertices_y[edge_vertex1];

        printf("%.4f %.4f %.4f %.4f\n", x0, y0, x1, y1);

        A += 0.5*(x0*y1-x1*y0);
        c_x_factor += (x0+x1)*(x0*y1-x1*y0);
        c_y_factor += (y0+y1)*(x0*y1-x1*y0);
      }

      mesh->cell_centers_x[cell_index] = (1.0/(6.0*A))*c_x_factor;
      mesh->cell_centers_y[cell_index] = (1.0/(6.0*A))*c_y_factor;


      printf("cell_centroids %d (%.4f %.4f)\n",
          cell_index,
          mesh->cell_centers_x[cell_index],
          mesh->cell_centers_y[cell_index]);
    }
  }
}

