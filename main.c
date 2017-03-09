#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "nodes_interface.h"
#include "unstructured_mesh.h"
#include "../profiler.h"
#include "../comms.h"
#include "../shared.h"
#include "../shared_data.h"
#include "../mesh.h"
#include "../params.h"

int main(int argc, char** argv) 
{
  if(argc < 4) {
    TERMINATE("Usage: ./nodes.exe <nx> <ny> <niters>\n");
  }

  Mesh mesh = {0};
  mesh.global_nx = atoi(argv[1]);
  mesh.global_ny = atoi(argv[2]);
  mesh.local_nx = atoi(argv[1]) + 2*PAD;
  mesh.local_ny = atoi(argv[2]) + 2*PAD;
  mesh.width = get_double_parameter("width", ARCH_ROOT_PARAMS);
  mesh.height = get_double_parameter("height", ARCH_ROOT_PARAMS);
  mesh.dt = get_double_parameter("max_dt", ARCH_ROOT_PARAMS);
  mesh.sim_end = get_double_parameter("sim_end", ARCH_ROOT_PARAMS);
  mesh.rank = MASTER;
  mesh.nranks = 1;
  mesh.niters = atoi(argv[3]);

  UnstructuredMesh unstructured_mesh;
  build_unstructured_mesh(
      &unstructured_mesh, mesh.global_nx, mesh.global_ny, mesh.local_nx, mesh.local_ny);

  initialise_mpi(argc, argv, &mesh.rank, &mesh.nranks);
  initialise_devices(mesh.rank);
  initialise_comms(&mesh);
  initialise_mesh_2d(&mesh);

  SharedData shared_data = {0};
  initialise_shared_data_2d(
      mesh.global_nx, mesh.global_ny, mesh.local_nx, mesh.local_ny, 
      mesh.x_off, mesh.y_off, &shared_data);

#if 0
  write_all_ranks_to_visit(
      mesh.global_nx+2*PAD, mesh.global_ny+2*PAD, mesh.local_nx, mesh.local_ny, mesh.x_off, 
      mesh.y_off, mesh.rank, mesh.nranks, mesh.neighbours, shared_data.x, "final_result", 0, 0.0);
#endif // if 0

  int tt = 0;
  double elapsed_sim_time = 0.0;
  double wallclock = 0.0;

  for(tt = 0; tt < mesh.niters; ++tt) {
    if(mesh.rank == MASTER) {
      printf("step %d\n", tt+1);
    }

    const double w0 = omp_get_wtime();

    int end_niters = 0;
    double end_error = 0.0;

    solve_unstructured_diffusion_2d(
        mesh.local_nx, mesh.local_ny, unstructured_mesh.nedges, 
        unstructured_mesh.vertices_x, unstructured_mesh.vertices_y, 
        unstructured_mesh.cells_vertices, 
        unstructured_mesh.nfaces, unstructured_mesh.density, 
        unstructured_mesh.cells_indirection1, unstructured_mesh.cells_indirection2, 
        unstructured_mesh.edges, unstructured_mesh.energy);

    wallclock += omp_get_wtime()-w0;

    if(mesh.rank == MASTER) {
      printf("finished on diffusion iteration %d with error %e\n", end_niters, end_error);
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
    printf("wallclock %.4f, elapsed simulation time %.4fs\n", wallclock, elapsed_sim_time);
  }

#if 0
  write_all_ranks_to_visit(
      mesh.global_nx+2*PAD, mesh.global_ny+2*PAD, mesh.local_nx, mesh.local_ny, mesh.x_off, 
      mesh.y_off, mesh.rank, mesh.nranks, mesh.neighbours, shared_data.p, "final_result", 0, elapsed_sim_time);
#endif // if 0

  finalise_shared_data(&shared_data);
  finalise_mesh(&mesh);
  finalise_comms();
}

