#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define NUMBER_OF_THREADS_PER_BLOCK 1024

__device__ double powCuda(double x, int) {
  return x * x;
}

int i4_ceiling(double x)
{
  int value = (int)x;
  if (value < x)
    value = value + 1;
  return value;
}

int i4_min(int i1, int i2)
{
  int value;
  if (i1 < i2)
    value = i1;
  else
    value = i2;
  return value;
}
__device__
double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (powCuda(x / a / a, 2) + powCuda(y / b / b, 2) + powCuda(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

__device__
double r8_uniform_01(int *seed)
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * (*seed - k * 127773) - k * 2836;

  if (*seed < 0)
  {
    *seed = *seed + 2147483647;
  }
  r = (double)(*seed) * 4.656612875E-10;

  return r;
}

void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}



__global__ void feymanKernel(int N, int ni, int nj, int nk, double* err, int* n_inside, double* wt) {

  int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

  const double a = 3.0;
  const double b = 2.0;
  const double c = 1.0;
  int i = blockIdx.x + 1;
  int j = blockIdx.y + 1;
  int k = blockIdx.z + 1;

  double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
  double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
  double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);

  double chk;
  double dx;
  double dy;
  double dz;
  const double h = 0.001;
  double stepsz;
  int seed = 123456789 + threadIdx.x;
  int steps;
  int trial;
  double us;
  double ut;
  double vh;
  double vs;
  double x1;
  double x2;
  double x3;
  double w;
  double w_exact;
  double we;
  steps = 0;

  const int dim = 3;
  stepsz = sqrt((double)dim * h);

  chk = powCuda(x / a, 2) + powCuda(y / b, 2) + powCuda(z / c, 2);
  
  if (1.0 < chk)
  {
    w_exact = 1.0;
    wt[blockId] = 1.0;
    return;
  }
  
  // Na nivou bloka, inkrementujemo n_inside
  if (threadIdx.x == 0) {
    atomicAdd(n_inside, 1);
    wt[blockId] = 0;
  }
  
  w_exact = exp(powCuda(x / a, 2) + powCuda(y / b, 2) + powCuda(z / c, 2) - 1.0);
  __syncthreads();
  double mywt = 0;
  for (trial = threadIdx.x; trial < N; trial += NUMBER_OF_THREADS_PER_BLOCK)
  {
    x1 = x;
    x2 = y;
    x3 = z;
    w = 1.0;
    chk = 0.0;
    while (chk < 1.0)
    {
      ut = r8_uniform_01(&seed);
      if (ut < 1.0 / 3.0)
      {
        us = r8_uniform_01(&seed) - 0.5;
        if (us < 0.0)
          dx = -stepsz;
        else
          dx = stepsz;
      }
      else
        dx = 0.0;

      ut = r8_uniform_01(&seed);
      if (ut < 1.0 / 3.0)
      {
        us = r8_uniform_01(&seed) - 0.5;
        if (us < 0.0)
          dy = -stepsz;
        else
          dy = stepsz;
      }
      else
        dy = 0.0;

      ut = r8_uniform_01(&seed);
      if (ut < 1.0 / 3.0)
      {
        us = r8_uniform_01(&seed) - 0.5;
        if (us < 0.0)
          dz = -stepsz;
        else
          dz = stepsz;
      }
      else
        dz = 0.0;

      vs = potential(a, b, c, x1, x2, x3);
      x1 = x1 + dx;
      x2 = x2 + dy;
      x3 = x3 + dz;

      steps++;

      vh = potential(a, b, c, x1, x2, x3);

      we = (1.0 - h * vs) * w;
      w = w - 0.5 * h * (vh * we + vs * w);

      chk = powCuda(x1 / a, 2) + powCuda(x2 / b, 2) + powCuda(x3 / c, 2);
    }
    // atomicAdd(&(wt[blockId]), w);  
    mywt += w;
  }

  atomicAdd(&(wt[blockId]), mywt); 
  __syncthreads();

  if (threadIdx.x == 0) {
    wt[blockId] = wt[blockId] / (double)(N);
  
    atomicAdd(err, powCuda(w_exact - wt[blockId], 2));
  }
}

// print na stdout upotrebiti u validaciji paralelnog resenja
int main(int arc, char **argv)
{
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  int dim = 3;
  double err;
  double h = 0.001;
  int n_inside;
  int ni;
  int nj;
  int nk;
  double stepsz;

  int N = atoi(argv[1]);
  timestamp();

  printf("A = %f\n", a);
  printf("B = %f\n", b);
  printf("C = %f\n", c);
  printf("N = %d\n", N);
  printf("H = %6.4f\n", h);

  stepsz = sqrt((double)dim * h);

  if (a == i4_min(i4_min(a, b), c))
  {
    ni = 6;
    nj = 1 + i4_ceiling(b / a) * (ni - 1);
    nk = 1 + i4_ceiling(c / a) * (ni - 1);
  }
  else if (b == i4_min(i4_min(a, b), c))
  {
    nj = 6;
    ni = 1 + i4_ceiling(a / b) * (nj - 1);
    nk = 1 + i4_ceiling(c / b) * (nj - 1);
  }
  else
  {
    nk = 6;
    ni = 1 + i4_ceiling(a / c) * (nk - 1);
    nj = 1 + i4_ceiling(b / c) * (nk - 1);
  }

  err = 0.0;
  n_inside = 0;

  dim3 dimGrid(ni, nj, nk);
  dim3 dimBlock(NUMBER_OF_THREADS_PER_BLOCK);

  int* device_n_inside;
  double* device_err;
  double* wt_device;

  int sharedMemSize = NUMBER_OF_THREADS_PER_BLOCK * sizeof(double);

  cudaMalloc((void**)&device_n_inside, sizeof(int));
  cudaMalloc((void**)&device_err, sizeof(double));
  cudaMalloc((void**)&wt_device, ni * nj * nk * sizeof(double));

  cudaMemcpy(device_n_inside, &n_inside, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(device_err, &err, sizeof(double), cudaMemcpyHostToDevice);

  // int N, double ni, double nj, double nk, double* err, int* n_inside

  feymanKernel<<< dimGrid, dimBlock >>>(N, ni, nj, nk, device_err, device_n_inside, wt_device); 
  cudaThreadSynchronize();

  cudaMemcpy(&n_inside, device_n_inside, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&err, device_err, sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(device_n_inside);
  cudaFree(device_err);
  cudaFree(wt_device);
  
  err = sqrt(err / (double)(n_inside));

  printf("\n\nRMS absolute error in solution = %e\n", err);
  timestamp();

  return 0;
}
