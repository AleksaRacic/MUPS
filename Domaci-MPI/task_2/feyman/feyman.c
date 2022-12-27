#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#define MASTER 0

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

double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

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

// print na stdout upotrebiti u validaciji paralelnog resenja
int main(int arc, char **argv)
{

  MPI_Init(&arc, &argv);

  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  double chk;
  int dim = 3;
  double dx;
  double dy;
  double dz;
  double err;
  double h = 0.001;
  int i;
  int j;
  int k;
  int n_inside;
  int ni;
  int nj;
  int nk;
  double stepsz;
  int seed = 123456789;
  int steps;
  int steps_ave;
  int trial;
  double us;
  double ut;
  double vh;
  double vs;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double w;
  double w_exact;
  double we;
  double wt;
  double z;
  int rank;
  int size;
  double master_wt;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  seed += rank;

  if(size > 4){
    MPI_Abort(MPI_COMM_WORLD, -1);
  }


  
  int N = atoi(argv[1]);

  double my_time = MPI_Wtime();

  if (rank == MASTER) {
    timestamp();

    printf("A = %f\n", a);
    printf("B = %f\n", b);
    printf("C = %f\n", c);
    printf("N = %d\n", N);
    printf("H = %6.4f\n", h);
  }

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

  int chunk_size = (N + size - 1) / size;
  int start = rank * chunk_size;
  int end = start + chunk_size > N ? N : start + chunk_size;

  double master_err;
  int master_n_inside;

  err = 0.0;
  n_inside = 0;

  //DEBUG
  int tot_trials = 0;
      
  for (i = 1; i <= ni; i++)
  {
    x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);

    for (j = 1; j <= nj; j++)
    {
      y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);

      for (k = 1; k <= nk; k++)
      {
        z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);

        chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

        if (1.0 < chk)
        {
          w_exact = 1.0;
          wt = 1.0;
          steps_ave = 0;
          // printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
          //        x, y, z, wt, w_exact, fabs(w_exact - wt), steps_ave);

          continue;
        }

        n_inside++;

        w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

        wt = 0.0;
        steps = 0;

        for (trial = start; trial < end; trial++)
        {
          //DEBUG
          tot_trials++;
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

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
          wt = wt + w;
        }
        MPI_Reduce(&wt, &master_wt, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        if(rank == MASTER){
          master_wt = master_wt / (double)(N);
          err = err + pow(w_exact - master_wt, 2);
        }
      }


    }
  }

  double maxtime, mintime, avgtime;

  my_time = MPI_Wtime() - my_time;

  MPI_Reduce(&my_time, &maxtime, 1, MPI_DOUBLE,MPI_MAX, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&my_time, &mintime, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
  MPI_Reduce(&my_time, &avgtime, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

  if (rank == MASTER) {
    err = sqrt(err / (double)(n_inside));

    printf("\n\nRMS absolute error in solution = %e\n", err);
    timestamp();

    avgtime /= size;

    FILE *fpt;
    fpt = fopen("feyman_MPI.csv", "a");
    fprintf(fpt,"%d, %d, %f, %f, %f, %f\n", size, N, err, maxtime, mintime, avgtime);
    fclose(fpt);
  }
  MPI_Finalize();

  return 0;
}
