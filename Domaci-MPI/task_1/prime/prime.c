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
  int k;
  int prime;
  int total;

  total = 0;
  
  int chunk_size = 0;
  int rank;
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 2 + 2 * rank; i <= n; i += 2 * size)
  {
    for(k = 0; k<2 && i + k <= n ; k++){
      prime = 1;
      for (j = 2; j < i + k; j++)
      {
        if (((i+k) % j) == 0)
        {
          prime = 0;
          break;
        }
    }
    total = total + prime;
    }
    
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

  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(size > 4){
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

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

  double maxtime, mintime, avgtime;
  double maxctime, minctime, avgctime;

  double my_time = MPI_Wtime();
  double my_ctime = cpu_time();

  test(n_lo, n_hi, n_factor);

  my_time = MPI_Wtime() - my_time;
  my_ctime = cpu_time() - my_ctime;

  //printf("\n--%d %f\n\n", rank, my_time);

  MPI_Reduce(&my_time, &maxtime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&my_time, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
  MPI_Reduce(&my_time, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

  MPI_Reduce(&my_ctime, &maxctime, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&my_ctime, &minctime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
  MPI_Reduce(&my_ctime, &avgctime, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

  if (rank == MASTER) {
    printf("\n");
    printf("PRIME_TEST\n");
    printf("  Normal end of execution.\n");
    printf("\n");
    timestamp();
    avgctime /= size;
    avgtime /= size;
    FILE *fpt;
    fpt = fopen("prime_parallel_MPI.csv", "a");
    fprintf(fpt,"%d, %f, %f, %f ,%f, %f, %f\n", size, maxtime, mintime, avgtime, maxctime, minctime, avgctime);
    fclose(fpt);
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
