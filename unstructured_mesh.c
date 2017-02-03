#include <stdio.h>
#include <stdlib.h>
#include "../mesh.h"
#include "../shared.h"
#include "../comms.h"
#include "unstructured_mesh.h"

/* Implementing a very manual approach to this in order to get something 
 * working and then it can be incorporated into the multi-package application
 * at some point */

// Build an unstructured mesh
void build_unstructured_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int local_nx, const int local_ny)
{
  allocate_data(&mesh->vertices_x, local_nx*local_ny);
  allocate_data(&mesh->vertices_y, local_nx*local_ny);
  allocate_data(&mesh->density, local_nx*local_ny);
  allocate_data(&mesh->energy, local_nx*local_ny);
  const int nedges = ((2*local_nx+1)*local_ny)-(local_nx+1);
  mesh->edges_vertices1 = (int*)malloc(sizeof(int)*nedges);
  mesh->edges_vertices2 = (int*)malloc(sizeof(int)*nedges);
  mesh->cells_indirection1 = (int*)malloc(sizeof(int)*nedges);
  mesh->cells_indirection2 = (int*)malloc(sizeof(int)*nedges);
  mesh->edges = (int*)malloc(sizeof(int)*local_nx*local_ny*NFACES);
  mesh->nfaces = (int*)malloc(sizeof(int)*local_nx*local_ny);

  /* Simply faking the unstructured mesh here, presumably it would be generated 
   * by some sort of parallel mesh generator */

  // Construct the list of vertices
  for(int ii = 0; ii < (local_ny+1); ++ii) {
    for(int jj = 0; jj < (local_nx+1); ++jj) {
      const int index = ii*(local_nx+1)+jj;
      mesh->vertices_x[index] = (double)(jj-PAD)*(WIDTH/(double)global_nx);
      mesh->vertices_y[index] = (double)(ii-PAD)*(HEIGHT/(double)global_ny);
    }
  }

  // Construct the list of edges
  // length of list (2*nx+1)*ny-nx
  int edge_index = 0;
  for(int ii = 0; ii < local_ny+1; ++ii) {
    // Add all of the edges along the x dimension for this row
    for(int jj = 0; jj < local_nx; ++jj) {
      mesh->edges_vertices1[edge_index] = (ii*(local_nx+1))+(jj);
      mesh->edges_vertices2[edge_index] = (ii*(local_nx+1))+(jj+1);
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
      mesh->edges_vertices1[edge_index] = ii*(local_nx+1)+jj;
      mesh->edges_vertices2[edge_index] = (ii+1)*(local_nx+1)+jj;
      mesh->cells_indirection1[edge_index] = 
        (jj > 0) ? (ii*local_nx)+(jj-1) : EDGE;
      mesh->cells_indirection2[edge_index] = 
        (jj < local_nx) ? (ii*local_nx)+(jj) : EDGE;
      edge_index++;
    }
  }

  // This artificial problem have four faces per cell
  for(int ii = 0; ii < local_nx*local_ny; ++ii) {
    mesh->nfaces[ii] = NFACES;
  }

  // Store the list of adjacent edges for each cell in order
  for(int ii = 0; ii < local_ny; ++ii) {
    for(int jj = 0; jj < local_nx; ++jj) {
      const int index = (ii*local_nx+jj)*NFACES;
      mesh->edges[index]   = ii*(2*local_nx+1)+jj;
      mesh->edges[index+1] = ii*(2*local_nx+1)+(jj+local_nx);
      mesh->edges[index+2] = (ii+1)*(2*local_nx+1)+jj;
      mesh->edges[index+3] = ii*(2*local_nx+1)+(jj+local_nx+1);
    }
  }

  printf("unstructured cartesian data initialised\n");
}

