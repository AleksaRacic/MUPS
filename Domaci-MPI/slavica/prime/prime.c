#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define MASTER 0

double cpu_time(void)
{
  double value;

  value = (double)clock() / (double)CLOCKS_PER_SEC;

  return value;
}

int prime_number(int n)
{
  int i;
  int j;
  int prime;
  int total;

  total = 0;
  
  int chunk_size = 0;
  int rank;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 2 + rank; i <= n; i += size)
  {
    prime = 1;
    for (j = 2; j < i; j++)
    {
      if ((i % j) == 0)
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  int master_total;
  MPI_Reduce(&total, &master_total, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);

  return master_total;
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
  MPI_Init(NULL, NULL);

  int n_factor;
  int n_hi;
  int n_lo;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == MASTER) {
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
  }
  
  MPI_Bcast(&n_lo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_hi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_factor, 1, MPI_INT, 0, MPI_COMM_WORLD);

  test(n_lo, n_hi, n_factor);

  if (rank == MASTER) {
    printf("\n");
    printf("PRIME_TEST\n");
    printf("  Normal end of execution.\n");
    printf("\n");
    timestamp();
  }

  MPI_Finalize();

  return 0;
}

void test(int n_lo, int n_hi, int n_factor)
{
  int i;
  int n;
  int primes;
  double ctime;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == MASTER) {
    printf("\n");
    printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
    printf("\n");
    printf("         N        Pi          Time\n");
    printf("\n");
  }

  n = n_lo;

  while (n <= n_hi)
  {
    if (rank == MASTER) {
      ctime = cpu_time();
    }

    primes = prime_number(n);

    if (rank == MASTER) {
      ctime = cpu_time() - ctime;

      printf("  %8d  %8d  %14f\n", n, primes, ctime);
    }

    n = n * n_factor;
  }

  return;
}
