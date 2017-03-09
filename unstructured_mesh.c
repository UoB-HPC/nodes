#include <stdio.h>
#include <stdlib.h>
#include "unstructured_mesh.h"
#include "../mesh.h"
#include "../shared.h"
#include "../comms.h"

/* Implementing a very manual approach to this in order to get something 
 * working and then it can be incorporated into the multi-package application
 * at some point */

// Build an unstructured mesh
void build_unstructured_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int local_nx, const int local_ny)
{
  allocate_data(&mesh->vertices_x, (local_nx+1)*(local_ny+1));
  allocate_data(&mesh->vertices_y, (local_nx+1)*(local_ny+1));
  allocate_data(&mesh->cell_centers_x, local_nx*local_ny);
  allocate_data(&mesh->cell_centers_y, local_nx*local_ny);
  allocate_data(&mesh->density, local_nx*local_ny);
  allocate_data(&mesh->energy, local_nx*local_ny);
  mesh->nedges = ((2*local_nx+1)*local_ny)-(local_nx+1);

  mesh->nfaces = (int*)malloc(sizeof(int)*local_nx*local_ny);
  mesh->cells_indirection1 = (int*)malloc(sizeof(int)*mesh->nedges);
  mesh->cells_indirection2 = (int*)malloc(sizeof(int)*mesh->nedges);
  mesh->edges = (int*)malloc(sizeof(int)*local_nx*local_ny*NFACES);

  // This artificial problem have four faces per cell
  for(int ii = 0; ii < local_nx*local_ny; ++ii) {
    mesh->nfaces[ii] = NFACES;
  }

  // Get the number of vertices attached to cells
  mesh->ncell_vertices = 0;
  for(int ii = 0; ii < local_ny; ++ii) {
    for(int jj = 0; jj < local_nx; ++jj) {
      const int cell_index = (ii*local_nx)+(jj);
      mesh->ncell_vertices += mesh->nfaces[cell_index];
    }
  }

  mesh->cells_vertices = (int*)malloc(sizeof(int)*mesh->ncell_vertices);

  /* Simply faking the unstructured mesh here, presumably it would be generated 
   * by some sort of parallel mesh generator */

  // Construct the list of vertices
  for(int ii = 0; ii < (local_ny+1); ++ii) {
    for(int jj = 0; jj < (local_nx+1); ++jj) {
      const int index = ii*(local_nx+1)+jj;
      mesh->vertices_x[index] = (double)(jj-PAD)*(mesh->width/(double)global_nx);
      mesh->vertices_y[index] = (double)(ii-PAD)*(mesh->height/(double)global_ny);
    }
  }

  int ncell_vertices = 0;
  for(int ii = 0; ii < local_ny; ++ii) {
    for(int jj = 0; jj < local_nx; ++jj) {
      const int cell_index = (ii*local_nx)+(jj);
      //for(int ii = 0; ii < mesh->nfaces[cell_index]; ++ii) {
        // The ordering here is clockwise and essential for later routines
        mesh->cells_vertices[ncell_vertices++] = (ii*(local_nx+1))+(jj);
        mesh->cells_vertices[ncell_vertices++] = (ii*(local_nx+1))+(jj+local_nx);
        mesh->cells_vertices[ncell_vertices++] = ((ii+1)*(local_nx+1))+(jj+local_nx+1);
        mesh->cells_vertices[ncell_vertices++] = ((ii+1)*(local_nx+1))+(jj);
      //}
    }
  }

  /* This operation can be done in a general way, by solving the equation of 
   * the lines created by the individual vertices crossing */

  // Find the (x,y) location of each of the cell centers
  int vertex_index = 0;
  for(int ii = 0; ii < local_ny; ++ii) {
    for(int jj = 0; jj < local_nx; ++jj) {
      const int cell_index = (ii*local_nx)+(jj);

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

  for(int ii = 0; ii < local_ny; ++ii) {
    for(int jj = 0; jj < local_nx; ++jj) {
      const int cell_index = (ii*local_nx)+(jj);
      printf("(%f, %f) ", 
          mesh->cell_centers_x[cell_index], mesh->cell_centers_y[cell_index]);
    }
    printf("\n");
  }

  printf("unstructured cartesian data initialised\n");
}


#if 0
int edge_index = 0;
for(int ii = 0; ii < local_ny+1; ++ii) {
  // Add all of the edges along the x dimension for this row
  for(int jj = 0; jj < local_nx; ++jj) {
    mesh->cells_indirection1[edge_index] = 
      (ii > 0) ? ((ii-1)*local_nx)+(jj) : EDGE;
    mesh->cells_indirection2[edge_index] = 
      (ii < local_ny) ? (ii*local_nx)+(jj) : EDGE;
    edge_index++;
  }
  // Don't fill in missing edges
  if(ii >= local_ny) {
    break;
  }
  // Add all of the edges along the y dimension for this row
  for(int jj = 0; jj < local_nx+1; ++jj) {
    mesh->cells_indirection1[edge_index] = 
      (jj > 0) ? (ii*local_nx)+(jj-1) : EDGE;
    mesh->cells_indirection2[edge_index] = 
      (jj < local_nx) ? (ii*local_nx)+(jj) : EDGE;
    edge_index++;
  }
}
#endif // if 0

