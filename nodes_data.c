#include <stdio.h>
#include <stdlib.h>
#include "nodes_data.h"
#include "nodes_interface.h"
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
  initialise_neighbour_list(nx, ny, nodes_data->neighbours_ii, nodes_data->neighbours_jj);
}

// Build an unstructured mesh
void build_unstructured_quad_mesh(
    UnstructuredMesh* mesh, const int global_nx, const int global_ny, 
    const int nx, const int ny)
{
  if(NFACES != 4) {
    TERMINATE("Only implemented for quads currently.");
  }

  mesh->nedges = 2*nx*ny+nx+ny;

  allocate_data(&mesh->vertices_x, (nx+1)*(ny+1));
  allocate_data(&mesh->vertices_y, (nx+1)*(ny+1));
  allocate_data(&mesh->density, nx*ny);
  allocate_data(&mesh->energy, nx*ny);

  // The area in each dimension of the face
  allocate_int_data(&mesh->face_area_x, mesh->nedges);
  allocate_int_data(&mesh->face_area_y, mesh->nedges);

  // Ordered by face_vertex_0 is closest to bottom left
  allocate_int_data(&mesh->face_vertex_0, mesh->nedges);
  allocate_int_data(&mesh->face_vertex_1, mesh->nedges);

  // Just setting all cells to have same number of faces
  const int nfaces = NFACES;

  // Construct the list of vertices contiguously
  for(int ii = 0; ii < (ny+1); ++ii) {
    for(int jj = 0; jj < (nx+1); ++jj) {
      const int index = (ii)*(nx+1)+(jj);
      mesh->vertices_x[index] = (double)((jj)-PAD)*(mesh->width/(double)global_nx);
      mesh->vertices_y[index] = (double)((ii)-PAD)*(mesh->height/(double)global_ny);
    }
  }

  // Calculate the vertices connecting each face, we step through from bottom 
  // left all the way to top right
  for(int ii = 0; ii < ny+1; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int face_index = ii*(2*nx+1)+jj;
      mesh->face_area_x[face_index] = 
        mesh->vertices_x[face_index+1]-mesh->vertices_x[face_index];
      mesh->face_area_y[face_indey] = 
        mesh->vertices_y[face_indey+1]-mesh->vertices_y[face_indey];
      mesh->face_vertex_0[face_index] = face_index;
      mesh->face_vertex_0[face_index] = face_index+1;
    }
    if(ii < ny) {
      for(int jj = 0; jj < nx+1; ++jj) {
        const int face_index = ii*(2*nx+1)+jj+nx;
        mesh->face_area_x[face_index] = 
          mesh->vertices_x[face_index+nx]-mesh->vertices_x[face_index];
        mesh->face_area_y[face_indey] = 
          mesh->vertices_y[face_indey+nx]-mesh->vertices_y[face_indey];
        mesh->face_vertex_0[face_index] = face_index;
        mesh->face_vertex_0[face_index] = face_index+nx;
      }
    }
  }

  allocate_int_data(&mesh->cells_faces, nx*ny);

  // Initialise cells connecting faces
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {

    }
  }

  allocate_data(&mesh->cell_centers_x, nx*ny);
  allocate_data(&mesh->cell_centers_y, nx*ny);

  // Find the (x,y) location of each of the cell centers
  int vertex_index = 0;
  for(int ii = 0; ii < ny; ++ii) {
    for(int jj = 0; jj < nx; ++jj) {
      const int cell_index = (ii*nx)+(jj);

      double A = 0.0;
      double c_x_factor = 0.0;
      double c_y_factor = 0.0;
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

#if 0
allocate_int_data(&mesh->nfaces, nx*ny);
for(int ii = 0; ii < nx*ny; ++ii) {
  mesh->nfaces[ii] = NFACES;
}
#endif // if 0

#if 0
// Get the number of vertices attached to cells
mesh->ncell_vertices = 0;
for(int ii = 0; ii < ny; ++ii) {
  for(int jj = 0; jj < nx; ++jj) {
    const int cell_index = (ii*nx)+(jj);
    mesh->ncell_vertices += mesh->nfaces[cell_index];
  }
}

mesh->cells_vertices = (int*)malloc(sizeof(int)*mesh->ncell_vertices);

/* This operation can be done in a general way, by solving the equation of 
 * the lines created by the individual vertices crossing */

for(int ii = 0; ii < ny; ++ii) {
  for(int jj = 0; jj < nx; ++jj) {
    const int cell_index = (ii*nx)+(jj);
    printf("(%f, %f) ", 
        mesh->cell_centers_x[cell_index], mesh->cell_centers_y[cell_index]);
  }
  printf("\n");
}

printf("unstructured cartesian data initialised\n");
#endif // if 0
