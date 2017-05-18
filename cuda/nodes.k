#include "../../cuda/reduction.k"

// https://inldigitallibrary.inl.gov/sti/3952796.pdf
// Take the average of the coefficients at the cells surrounding 
// each face
__global__ void calc_s_x(
    const int nx, const int ny, const double dt, const double heat_capacity,
    const double conductivity, double* s_x, const double* rho, const double* edgedx) 
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx+1));
  const int ii = (gid / (nx+1));
  if(ii < PAD || jj < PAD || jj >= (nx+1)-PAD || ii >= ny-PAD) return;
  s_x[(ii)*(nx+1)+(jj)] = (dt*conductivity*(rho[(ii)*nx+(jj)]+rho[(ii)*nx+(jj-1)]))/
    (2.0*rho[(ii)*nx+(jj)]*rho[(ii)*nx+(jj-1)]*edgedx[jj]*edgedx[jj]*heat_capacity);
}

__global__ void calc_s_y(
    const int nx, const int ny, const double dt, const double heat_capacity,
    const double conductivity, double* s_y, const double* rho, const double* edgedy) 
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx));
  const int ii = (gid / (nx));
  if(ii < PAD || jj < PAD || jj >= nx-PAD || ii >= (ny+1)-PAD) return;
  s_y[(ii)*nx+(jj)] = (dt*conductivity*(rho[(ii)*nx+(jj)]+rho[(ii-1)*nx+(jj)]))/
    (2.0*rho[(ii)*nx+(jj)]*rho[(ii-1)*nx+(jj)]*edgedy[ii]*edgedy[ii]*heat_capacity);
}

__global__ void calc_initial_r2(
    const int nx, const int ny, double* s_x, double* s_y, const double* x, 
    double* p, double* r, double* reduce_array, const int nneighbours, 
    const int* neighbours_ii, const int* neighbours_jj) 
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx));
  const int ii = (gid / (nx));
  const int tid = threadIdx.x;

  __shared__ double sdata[NTHREADS];
  sdata[tid] = 0.0;

  if(ii >= PAD && jj >= PAD && jj < nx-PAD && ii < ny-PAD) {
    const int neighbour_index = (ii)*nx*nneighbours+(jj)*nneighbours;
    const int nii = neighbours_ii[neighbour_index+NORTH_STENCIL];
    const int eii = neighbours_ii[neighbour_index+EAST_STENCIL];
    const int sii = neighbours_ii[neighbour_index+SOUTH_STENCIL];
    const int wii = neighbours_ii[neighbour_index+WEST_STENCIL];
    const int njj = neighbours_jj[neighbour_index+NORTH_STENCIL];
    const int ejj = neighbours_jj[neighbour_index+EAST_STENCIL];
    const int sjj = neighbours_jj[neighbour_index+SOUTH_STENCIL];
    const int wjj = neighbours_jj[neighbour_index+WEST_STENCIL];

    r[(ii)*nx+(jj)] = x[(ii)*nx+(jj)] -
      ((s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
      s_x[(eii)*(nx+1)+(ejj)]+s_y[(nii)*nx+(njj)])*x[(ii)*nx+(jj)]
       - s_y[(ii)*nx+(jj)]*x[(sii)*nx+(sjj)]
       - s_x[(ii)*(nx+1)+(jj)]*x[(wii)*nx+(wjj)] 
       - s_x[(eii)*(nx+1)+(ejj)]*x[(eii)*nx+(ejj)]
       - s_y[(nii)*nx+(jj)]*x[(nii)*nx+(jj)]);
    p[(ii)*nx+(jj)] = r[(ii)*nx+(jj)];
    sdata[tid] = r[(ii)*nx+(jj)]*r[(ii)*nx+(jj)];
  }
  __syncthreads();

  sum_reduce_in_shared<double, NTHREADS>(threadIdx.x, sdata);
  if(tid == 0) 
    reduce_array[blockIdx.x] = sdata[0];
}

__global__ void calc_pAp(
    const int nx, const int ny, const double* s_x, const double* s_y, double* p, 
    double* Ap, double* reduce_array, const int nneighbours, 
    const int* neighbours_ii, const int* neighbours_jj)
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx));
  const int ii = (gid / (nx));
  const int tid = threadIdx.x;

  __shared__ double sdata[NTHREADS];
  sdata[tid] = 0.0;

  // You don't need to use a matrix as the band matrix is fully predictable
  // from the 5pt stencil
  if(ii >= PAD & jj >= PAD && jj < nx-PAD && ii < ny-PAD) {
    const int neighbour_index = (ii)*nx*nneighbours+(jj)*nneighbours;
    const int nii = neighbours_ii[neighbour_index+NORTH_STENCIL];
    const int eii = neighbours_ii[neighbour_index+EAST_STENCIL];
    const int sii = neighbours_ii[neighbour_index+SOUTH_STENCIL];
    const int wii = neighbours_ii[neighbour_index+WEST_STENCIL];
    const int njj = neighbours_jj[neighbour_index+NORTH_STENCIL];
    const int ejj = neighbours_jj[neighbour_index+EAST_STENCIL];
    const int sjj = neighbours_jj[neighbour_index+SOUTH_STENCIL];
    const int wjj = neighbours_jj[neighbour_index+WEST_STENCIL];

    Ap[(ii)*nx+(jj)] = 
      (s_y[(ii)*nx+(jj)]+s_x[(ii)*(nx+1)+(jj)]+1.0+
       s_x[(eii)*(nx+1)+(ejj)]+s_y[(nii)*nx+(njj)])*p[(ii)*nx+(jj)]
      - s_y[(ii)*nx+(jj)]*p[(sii)*nx+(sjj)]
      - s_x[(ii)*(nx+1)+(jj)]*p[(wii)*nx+(wjj)] 
      - s_x[(eii)*(nx+1)+(ejj)]*p[(eii)*nx+(ejj)]
      - s_y[(nii)*nx+(njj)]*p[(nii)*nx+(njj)];
    sdata[threadIdx.x] = p[(ii)*nx+(jj)]*Ap[(ii)*nx+(jj)];
  }
  __syncthreads();

  sum_reduce_in_shared<double, NTHREADS>(threadIdx.x, sdata);
  if (tid == 0) 
    reduce_array[blockIdx.x] = sdata[0];
}

__global__ void calc_new_r2(
    const int nx, const int ny, double alpha, 
    double* x, double* p, double* r, double* Ap, double* reduce_array)
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx));
  const int ii = (gid / (nx));
  const int tid = threadIdx.x;

  __shared__ double sdata[NTHREADS];
  sdata[tid] = 0.0;

  if(ii >= PAD & jj >= PAD && jj < nx-PAD && ii < ny-PAD) {
    x[(ii)*nx+(jj)] += alpha*p[(ii)*nx+(jj)];
    r[(ii)*nx+(jj)] -= alpha*Ap[(ii)*nx+(jj)];
    sdata[threadIdx.x] = r[(ii)*nx+(jj)]*r[(ii)*nx+(jj)];
  }
  __syncthreads();

  sum_reduce_in_shared<double, NTHREADS>(threadIdx.x, sdata);
  if (tid == 0) 
    reduce_array[blockIdx.x] = sdata[0];
}

__global__ void update_p(
    const int nx, const int ny, const double beta, const double* r, double* p)
{
  const int gid = threadIdx.x+blockIdx.x*blockDim.x;
  const int jj = (gid % (nx));
  const int ii = (gid / (nx));
  if(ii < PAD || jj < PAD || jj >= nx-PAD || ii >= ny-PAD) return;

  p[(ii)*nx+(jj)] = r[(ii)*nx+(jj)] + beta*p[(ii)*nx+(jj)];
}

