
/*
 *  Compute forces and accumulate the virial and the potential
 */
#include <mpi.h>

#define MASTER 0

extern double epot, vir;


void forces(int npart, double x[], double f[], double side, double rcoff)
{
  int i, j;
  double sideh, rcoffs;
  double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
  double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
  double forcex, forcey, forcez;

  vir = 0.0;
  epot = 0.0;
  sideh = 0.5 * side;
  rcoffs = rcoff * rcoff;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == MASTER) {
    int proc_size;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    int size = proc_size - 1;
    
    for (i = 1;i<proc_size;i++) {
      int chunk_size = (npart * 3 + size - 1) / size;
      MPI_Send( &chunk_size, 1, MPI_INT, i, i, MPI_COMM_WORLD );
    }

    MPI_Status status;
    double tmp_f[npart * 3];
    double tmp_vir;
    double tmp_epot;
    for (i = 1; i < proc_size; i++) {
      MPI_Recv(tmp_f, npart * 3, MPI_DOUBLE, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, &status );
      MPI_Recv(&tmp_epot, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 101, MPI_COMM_WORLD, &status );
      MPI_Recv(&tmp_vir, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 102, MPI_COMM_WORLD, &status );

      for (j = 0; j < npart * 3; j++) {
        f[j] += tmp_f[j];
      }

      vir += tmp_vir;
      epot += tmp_epot;
    }

  } else {
    MPI_Status status;
    int chunk_size;
    MPI_Recv( &chunk_size, 1, MPI_INT, MASTER, rank, MPI_COMM_WORLD, &status );

    int start = (rank - 1) * chunk_size;
    int end = start + chunk_size > npart * 3 ? npart * 3 : start + chunk_size;

    for (i = start; i < end; i += 3)
    {
      xi = x[i];
      yi = x[i + 1];
      zi = x[i + 2];
      fxi = 0.0;
      fyi = 0.0;
      fzi = 0.0;

      for (j = i + 3; j < npart * 3; j += 3)
      {
        xx = xi - x[j];
        yy = yi - x[j + 1];
        zz = zi - x[j + 2];
        if (xx < -sideh)
          xx += side;
        if (xx > sideh)
          xx -= side;
        if (yy < -sideh)
          yy += side;
        if (yy > sideh)
          yy -= side;
        if (zz < -sideh)
          zz += side;
        if (zz > sideh)
          zz -= side;
        rd = xx * xx + yy * yy + zz * zz;

        if (rd <= rcoffs)
        {
          rrd = 1.0 / rd;
          rrd2 = rrd * rrd;
          rrd3 = rrd2 * rrd;
          rrd4 = rrd2 * rrd2;
          rrd6 = rrd2 * rrd4;
          rrd7 = rrd6 * rrd;
          epot += (rrd6 - rrd3);
          r148 = rrd7 - 0.5 * rrd4;
          vir -= rd * r148;
          forcex = xx * r148;
          fxi += forcex;
          f[j] -= forcex;
          forcey = yy * r148;
          fyi += forcey;
          f[j + 1] -= forcey;
          forcez = zz * r148;
          fzi += forcez;
          f[j + 2] -= forcez;
        }
      }
      f[i] += fxi;
      f[i + 1] += fyi;
      f[i + 2] += fzi;
    }

    MPI_Send(f, npart * 3, MPI_DOUBLE, MASTER, 100, MPI_COMM_WORLD);
    MPI_Send(&epot, 1, MPI_DOUBLE, MASTER, 101, MPI_COMM_WORLD);
    MPI_Send(&vir, 1, MPI_DOUBLE, MASTER, 102, MPI_COMM_WORLD);
  }
}
