#include <stdio.h>
#include <stdlib.h>
#include "nodes_data.h"
#include "../params.h"
#include "../mesh.h"
#include "../shared.h"
#include "../comms.h"

/* Implementing a very manual approach to this in order to get something 
 * working and then it can be incorporated into the multi-package application
 * at some point */

void initialise_nodes_data(
    NodesData* nodes_data, const int nx, const int ny, const int nneighbours, 
    const char* nodes_params)
{
  nodes_data->conductivity = get_double_parameter("conductivity", nodes_params);
  nodes_data->heat_capacity = get_double_parameter("heat_capacity", nodes_params);
  nodes_data->nneighbours = nneighbours;

  // Describe the neighbours for all real cells
  allocate_int_data(&nodes_data->neighbours_ii, nneighbours*nx*ny);
  allocate_int_data(&nodes_data->neighbours_jj, nneighbours*nx*ny);

  int* h_neighbours_ii;
  int* h_neighbours_jj;
  allocate_host_int_data(&h_neighbours_ii, nneighbours*nx*ny);
  allocate_host_int_data(&h_neighbours_jj, nneighbours*nx*ny);

  for(int ii = PAD; ii < ny-PAD; ++ii) {
    for(int jj = PAD; jj < nx-PAD; ++jj) {
      // Manual 5pt stencil
      h_neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+NORTH_STENCIL] = (ii+1);
      h_neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+EAST_STENCIL] = (ii);
      h_neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+SOUTH_STENCIL] = (ii-1);
      h_neighbours_ii[(ii)*nx*nneighbours+(jj)*nneighbours+WEST_STENCIL] = (ii);
      h_neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+NORTH_STENCIL] = (jj);
      h_neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+EAST_STENCIL] = (jj+1);
      h_neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+SOUTH_STENCIL] = (jj);
      h_neighbours_jj[(ii)*nx*nneighbours+(jj)*nneighbours+WEST_STENCIL] = (jj-1);
    }
  }

  // TODO: FIX THE STUPID MEMORY LEAK HERE
  copy_int_buffer(nneighbours*nx*ny, &h_neighbours_ii, &nodes_data->neighbours_ii, SEND);
  copy_int_buffer(nneighbours*nx*ny, &h_neighbours_jj, &nodes_data->neighbours_jj, SEND);
}

#if 0
// Build an unstructured mesh
void build_unstructured_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int nx, const int ny)
{
  allocate_data(&mesh->vertices_x, (nx+1)*(ny+1));
  allocate_data(&mesh->vertices_y, (nx+1)*(ny+1));
  allocate_data(&mesh->cell_centers_x, nx*ny);
  allocate_data(&mesh->cell_centers_y, nx*ny);
  allocate_data(&mesh->density, nx*ny);
  allocate_data(&mesh->energy, nx*ny);
  mesh->nedges = ((2*nx+1)*ny)-(nx+1);

  mesh->nfaces = (int*)malloc(sizeof(int)*nx*ny);
  mesh->cells_indirection1 = (int*)malloc(sizeof(int)*mesh->nedges);
  mesh->cells_indirection2 = (int*)malloc(sizeof(int)*mesh->nedges);
  mesh->edges = (int*)malloc(sizeof(int)*nx*ny*NFACES);

  // This artificial problem have four faces per cell
  for(int ii = 0; ii < nx*ny; ++ii) {
    mesh->nfaces[ii] = NFACES;
  }

  // Get the number of vertices attached to cells
  mesh->ncell_vertices = 0;
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii*nx)+(jj);
      mesh->ncell_vertices += mesh->nfaces[cell_index];
    }
  }

  mesh->cells_vertices = (int*)malloc(sizeof(int)*mesh->ncell_vertices);

  /* Simply faking the unstructured mesh here, presumably it would be generated 
   * by some sort of parallel mesh generator */

  // Construct the list of vertices
  for(int ii = 0; ii < (ny+1); ++ii) {
    for(int jj = 0; jj < (nx+1); ++jj) {
      const int index = ii*(nx+1)+jj;
      mesh->vertices_x[index] = (double)(jj-PAD)*(mesh->width/(double)global_nx);
      mesh->vertices_y[index] = (double)(ii-PAD)*(mesh->height/(double)global_ny);
    }
  }

  int ncell_vertices = 0;
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii*nx)+(jj);
      //for(int ii = 0; ii < mesh->nfaces[cell_index]; ++ii) {
      // The ordering here is clockwise and essential for later routines
      mesh->cells_vertices[ncell_vertices++] = (ii*(nx+1))+(jj);
      mesh->cells_vertices[ncell_vertices++] = (ii*(nx+1))+(jj+nx);
      mesh->cells_vertices[ncell_vertices++] = ((ii+1)*(nx+1))+(jj+nx+1);
      mesh->cells_vertices[ncell_vertices++] = ((ii+1)*(nx+1))+(jj);
      //}
    }
  }

  /* This operation can be done in a general way, by solving the equation of 
   * the lines created by the individual vertices crossing */

  // Find the (x,y) location of each of the cell centers
  int vertex_index = 0;
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii*nx)+(jj);

      double A = 0.0;
      double c_x_factor = 0.0;
      double c_y_factor = 0.0;
      const int nfaces = mesh->nfaces[cell_index];
      for(int kk = 0; kk < nfaces; ++kk) {
        const int vertex_index0 = 
          mesh->cells_vertices[vertex_index];
        const int vertex_index1 = 
          mesh->cells_vertices[(vertex_index+1)%nfaces];
        const double x0 = mesh->vertices_x[vertex_index0];
        const double y0 = mesh->vertices_y[vertex_index0];
        const double x1 = mesh->vertices_x[vertex_index1];
        const double y1 = mesh->vertices_y[vertex_index1];

        A += (x0*y1-x1*y0);
        c_x_factor += (x0+x1)*(x0*y1-x1*y0);
        c_y_factor += (y0+y1)*(x0*y1-x1*y0);

        // Nasty thread unsafe
        vertex_index++;
      }

      mesh->cell_centers_x[cell_index] = (1.0/(6.0*0.5*A))*c_x_factor;
      mesh->cell_centers_y[cell_index] = (1.0/(6.0*0.5*A))*c_y_factor;
    }
  }

  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii*nx)+(jj);
      printf("(%f, %f) ", 
          mesh->cell_centers_x[cell_index], mesh->cell_centers_y[cell_index]);
    }
    printf("\n");
  }

  printf("unstructured cartesian data initialised\n");
}


#if 0
int edge_index = 0;
for(int ii = 0; ii < ny+1; ++ii) {
  // Add all of the edges along the x dimension for this row
  for(int jj = 0; jj < nx; ++jj) {
    mesh->cells_indirection1[edge_index] = 
      (ii > 0) ? ((ii-1)*nx)+(jj) : EDGE;
    mesh->cells_indirection2[edge_index] = 
      (ii < ny) ? (ii*nx)+(jj) : EDGE;
    edge_index++;
  }
  // Don't fill in missing edges
  if(ii >= ny) {
    break;
  }
  // Add all of the edges along the y dimension for this row
  for(int jj = 0; jj < nx+1; ++jj) {
    mesh->cells_indirection1[edge_index] = 
      (jj > 0) ? (ii*nx)+(jj-1) : EDGE;
    mesh->cells_indirection2[edge_index] = 
      (jj < nx) ? (ii*nx)+(jj) : EDGE;
    edge_index++;
  }
}
#endif // if 0
#endif // if 0

