#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include<sys/time.h>

double wall_time(void) {
    struct timeval tv;

    gettimeofday(&tv,NULL);

    return (double)tv.tv_usec / 1000000 + (double)tv.tv_sec;
}

double cpu_time(void)
{
  double value;

  value = (double)clock() / (double)CLOCKS_PER_SEC;

  return value;
}

__global__ void primeNumberKernel(int n, int *total)
{
  int j;
  int number = 2 + blockIdx.x * blockDim.x + threadIdx.x;

  if (number > n)
    return;

  for (j = 2; j < number; j ++)
  {
    if ((number % j) == 0) return;
  }

  atomicAdd(total, 1);
}

#define NUMBER_OF_THREADS_PER_BLOCK 1024

int prime_number(int n)
{
  int total;
  total = 0;

  int numBlock = (n + NUMBER_OF_THREADS_PER_BLOCK - 1) / NUMBER_OF_THREADS_PER_BLOCK;

  dim3 dimGrid(numBlock);
  dim3 dimBlock(NUMBER_OF_THREADS_PER_BLOCK);

  int *total_device;
  cudaMalloc((void **)&total_device, sizeof(int));
  cudaMemcpy(total_device, &total, sizeof(int), cudaMemcpyHostToDevice);

  primeNumberKernel<<<dimGrid, dimBlock>>>(n, total_device);

  cudaDeviceSynchronize();

  cudaMemcpy(&total, total_device, sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(total_device);

  return total;
}

void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}

void test(int n_lo, int n_hi, int n_factor);

int main(int argc, char *argv[])
{
  int n_factor;
  int n_hi;
  int n_lo;

  double w_time;

  timestamp();
  printf("\n");
  printf("PRIME TEST\n");

  if (argc != 4)
  {
    n_lo = 1;
    n_hi = 131072;
    n_factor = 2;
  }
  else
  {
    n_lo = atoi(argv[1]);
    n_hi = atoi(argv[2]);
    n_factor = atoi(argv[3]);
  }

  w_time = wall_time();
  test(n_lo, n_hi, n_factor);
  w_time = wall_time() - w_time;
  
  printf("\n");
  printf("PRIME_TEST\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  FILE *fpt;
  fpt = fopen("task_1.csv", "a");
  fprintf(fpt,"CUDA, %f\n", w_time);
  fclose(fpt);

  return 0;
}

void test(int n_lo, int n_hi, int n_factor)
{
  int n;
  int primes;
  double ctime;

  printf("\n");
  printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
  printf("\n");
  printf("         N        Pi          Time\n");
  printf("\n");

  n = n_lo;

  while (n <= n_hi)
  {
    ctime = cpu_time();

    primes = prime_number(n);

    ctime = cpu_time() - ctime;

    printf("  %8d  %8d  %14f\n", n, primes, ctime);
    n = n * n_factor;
  }

  return;
}
