#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "nodes_interface.h"
#include "nodes_data.h"
#include "../profiler.h"
#include "../comms.h"
#include "../shared.h"
#include "../shared_data.h"
#include "../mesh.h"
#include "../params.h"

int main(int argc, char** argv) 
{
  if(argc < 2) {
    TERMINATE("Usage: ./nodes.exe <parameter_filename>\n");
  }

  Mesh mesh = {0};
  const char* nodes_params = argv[1];
  mesh.global_nx = get_int_parameter("nx", nodes_params);
  mesh.global_ny = get_int_parameter("ny", nodes_params);
  mesh.local_nx = mesh.global_nx + 2*PAD;
  mesh.local_ny = mesh.global_ny + 2*PAD;
  mesh.width = get_double_parameter("width", ARCH_ROOT_PARAMS);
  mesh.height = get_double_parameter("height", ARCH_ROOT_PARAMS);
  mesh.sim_end = get_double_parameter("sim_end", ARCH_ROOT_PARAMS);
  mesh.dt = get_double_parameter("max_dt", ARCH_ROOT_PARAMS);
  mesh.rank = MASTER;
  mesh.nranks = 1;
  mesh.niters = get_int_parameter("iterations", nodes_params);
  const int max_inners = get_int_parameter("max_inners", nodes_params);
  const int visit_dump = get_int_parameter("visit_dump", nodes_params);

  initialise_mpi(argc, argv, &mesh.rank, &mesh.nranks);
  initialise_devices(mesh.rank);
  initialise_comms(&mesh);
  initialise_mesh_2d(&mesh);

  UnstructuredMesh unstructured_mesh;
  initialise_unstructured_quad_mesh_2d(&unstructured_mesh, &mesh);

  NodesData nodes_data = {0};
  initialise_nodes_data(
      mesh.local_nx, mesh.local_ny, &nodes_data, nodes_params);

  SharedData shared_data = {0};
  initialise_shared_data_2d(
      mesh.global_nx, mesh.global_ny, mesh.local_nx, mesh.local_ny, mesh.x_off, 
      mesh.y_off, mesh.width, mesh.height, nodes_params, mesh.edgex, 
      mesh.edgey, &shared_data);

  handle_boundary_2d(
      mesh.local_nx, mesh.local_ny, &mesh, shared_data.rho, NO_INVERT, PACK);
  handle_boundary_2d(
      mesh.local_nx, mesh.local_ny, &mesh, shared_data.e, NO_INVERT, PACK);
  handle_boundary_2d(
      mesh.local_nx, mesh.local_ny, &mesh, shared_data.x, NO_INVERT, PACK);

  if(visit_dump) {
    write_all_ranks_to_visit(
        mesh.global_nx+2*PAD, mesh.global_ny+2*PAD, mesh.local_nx, mesh.local_ny, 
        mesh.x_off, mesh.y_off, mesh.rank, mesh.nranks, mesh.neighbours, 
        shared_data.x, "final_result", 0, 0.0);
  }

  int tt = 0;
  double elapsed_sim_time = 0.0;
  double wallclock = 0.0;

  for(tt = 0; tt < mesh.niters; ++tt) {
    if(mesh.rank == MASTER) {
      printf("step %d\n", tt+1);
    }

    double w0 = omp_get_wtime();

    int end_niters = 0;
    double end_error = 0.0;
    solve_unstructured_diffusion_2d(
        mesh.local_nx, mesh.local_ny, &mesh, &unstructured_mesh, max_inners, mesh.dt, 
        nodes_data.heat_capacity, nodes_data.conductivity, shared_data.x, 
        nodes_data.b, shared_data.r, shared_data.p, shared_data.rho, shared_data.Ap, 
        &end_niters, &end_error, shared_data.reduce_array0);

    wallclock += omp_get_wtime()-w0;

    if(mesh.rank == MASTER) {
      printf("finished on diffusion iteration %d with error %e\n", 
          end_niters, end_error);
    }

    elapsed_sim_time += mesh.dt;
    if(elapsed_sim_time >= mesh.sim_end) {
      if(mesh.rank == MASTER) {
        printf("reached end of simulation time\n");
      }
      break;
    }
  }

  if(mesh.rank == MASTER) {
    PRINT_PROFILING_RESULTS(&compute_profile);
    printf("wallclock %.4fs, elapsed simulation time %.4fs\n", 
        wallclock, elapsed_sim_time);
  }

  if(visit_dump) {
    write_all_ranks_to_visit(
        mesh.global_nx+2*PAD, mesh.global_ny+2*PAD, mesh.local_nx, mesh.local_ny, 
        mesh.x_off, mesh.y_off, mesh.rank, mesh.nranks, mesh.neighbours, 
        shared_data.x, "final_result", 1, elapsed_sim_time);
  }




  write_curvilinear_data_to_visit(
      mesh.local_nx, mesh.local_ny, 0, unstructured_mesh.vertices_x, 
      unstructured_mesh.vertices_y, shared_data.x);




  finalise_shared_data(&shared_data);
  finalise_mesh(&mesh);
  finalise_comms();
}

