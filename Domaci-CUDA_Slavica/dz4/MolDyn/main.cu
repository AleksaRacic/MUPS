#include <stdio.h>
#include <math.h>
#include <time.h>
// #include <omp.h>

#define NUMBER_OF_THREADS_PER_BLOCK 1024

#define mm 15
#define npart 4 * mm *mm *mm
/*
 *  Function declarations
 */

void dfill(int, double, double[], int);

void domove(int, double[], double[], double[], double);

void dscal(int, double, double[], int);

void fcc(double[], int, int, double);

__global__ void forces(int, double*, double*, double, double, double*, double*);

double
mkekin(int, double[], double[], double, double);

void mxwell(double[], int, double, double);

void prnout(int, double, double, double, double, double, double, int, double);

double
velavg(int, double[], double, double);

double
secnds(void);

/*
 *  Variable declarations
 */

double epot;
double vir;
double count;

/*
 *  Main program : Molecular Dynamics simulation.
 */
int main()
{
  int move;
  double x[npart * 3], vh[npart * 3], f[npart * 3];
  double ekin;
  double vel;
  double sc;
  double start, time;

  /*
   *  Parameter definitions
   */

  double den = 0.83134;
  double side = pow((double)npart / den, 0.3333333);
  double tref = 0.722;
  double rcoff = (double)mm / 4.0;
  double h = 0.064;
  int irep = 10;
  int istop = 20;
  int iprint = 5;
  int movemx = 20;

  double a = side / (double)mm;
  double hsq = h * h;
  double hsq2 = hsq * 0.5;
  double tscale = 16.0 / ((double)npart - 1.0);
  double vaver = 1.13 * sqrt(tref / 24.0);

  /*
   *  Initial output
   */

  printf(" Molecular Dynamics Simulation example program\n");
  printf(" ---------------------------------------------\n");
  printf(" number of particles is ............ %6d\n", npart);
  printf(" side length of the box is ......... %13.6f\n", side);
  printf(" cut off is ........................ %13.6f\n", rcoff);
  printf(" reduced temperature is ............ %13.6f\n", tref);
  printf(" basic timestep is ................. %13.6f\n", h);
  printf(" temperature scale interval ........ %6d\n", irep);
  printf(" stop scaling at move .............. %6d\n", istop);
  printf(" print interval .................... %6d\n", iprint);
  printf(" total no. of steps ................ %6d\n", movemx);

  /*
   *  Generate fcc lattice for atoms inside box
   */
  fcc(x, npart, mm, a);
  /*
   *  Initialise velocities and forces (which are zero in fcc positions)
   */
  mxwell(vh, 3 * npart, h, tref);
  dfill(3 * npart, 0.0, f, 1);
  /*
   *  Start of md
   */
  printf("\n    i       ke         pe            e         temp   "
         "   pres      vel      rp\n  -----  ----------  ----------"
         "  ----------  --------  --------  --------  ----\n");

  start = secnds();

  for (move = 1; move <= movemx; move++)
  {

    /*
     *  Move the particles and partially update velocities
     */
    domove(3 * npart, x, vh, f, side);

    /*
     *  Compute forces in the new positions and accumulate the virial
     *  and potential energy.
     */

    int blockSize = (npart - 1 + NUMBER_OF_THREADS_PER_BLOCK)/NUMBER_OF_THREADS_PER_BLOCK;

    dim3 dimGrid(blockSize);
    dim3 dimBlock(NUMBER_OF_THREADS_PER_BLOCK);

    int sharedMemSize = NUMBER_OF_THREADS_PER_BLOCK * npart * sizeof(double);
    double* f_device;
    double* x_device;
    double* epot_device;
    double* vir_device;
    vir = 0.0;
    epot = 0.0;
    cudaMalloc((void**)&f_device,npart * 3 * sizeof(double));
    cudaMalloc((void**)&x_device,npart * 3 * sizeof(double));
    cudaMalloc((void**)&epot_device,sizeof(double));
    cudaMalloc((void**)&vir_device,sizeof(double));

    double npart_res = 4 * mm *mm *mm;

    cudaMemcpy(f_device, f, npart * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(x_device, x, npart * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(epot_device, &epot, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(vir_device, &vir, sizeof(double), cudaMemcpyHostToDevice);

      //(int *npart, double *x, double *f, double *side, double *rcoff, double *vir, double *epot)
    forces<<< dimGrid, dimBlock >>>(npart, x_device, f_device, side, rcoff, vir_device, epot_device);

    cudaThreadSynchronize();

    cudaMemcpy(f, f_device, npart * 3 * sizeof(double), cudaMemcpyDeviceToHost);

    cudaMemcpy(x, x_device, npart * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&epot, epot_device, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&vir, vir_device, sizeof(double), cudaMemcpyDeviceToHost);

    // printf("\n%lf    %lf\n", epot, vir);

    cudaFree(f_device);
    cudaFree(x_device);
    cudaFree(epot_device);
    cudaFree(vir_device);
    

    /*
     *  Scale forces, complete update of velocities and compute k.e.
     */
    ekin = mkekin(npart, f, vh, hsq2, hsq);

    /*
     *  Average the velocity and temperature scale if desired
     */
    vel = velavg(npart, vh, vaver, h);
    if (move < istop && fmod(move, irep) == 0)
    {
      sc = sqrt(tref / (tscale * ekin));
      dscal(3 * npart, sc, vh, 1);
      ekin = tref / tscale;
    }

    /*
     *  Sum to get full potential energy and virial
     */
    if (fmod(move, iprint) == 0)
      prnout(move, ekin, epot, tscale, vir, vel, count, npart, den);
  }

  time = secnds() - start;

  printf("Time =  %f\n", (float)time);
}

time_t starttime = 0;

double secnds()
{

  return time(0);//omp_get_wtime();
}
