
/*
 *  Compute forces and accumulate the virial and the potential
 */


__global__ void forces(int npart, double *x, double *f, double side, double rcoff, double *vir, double *epot)
{
  int i, j;

  i = (blockIdx.x * blockDim.x + threadIdx.x) * 3;

  double sideh, rcoffs;
  double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
  double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
  double forcex, forcey, forcez;

  // vir = 0.0;
  // epot = 0.0;
  sideh = 0.5 * side;
  rcoffs = rcoff * rcoff;

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
      // epot += (rrd6 - rrd3);
      atomicAdd(epot, (rrd6 - rrd3));
      r148 = rrd7 - 0.5 * rrd4;
      // vir -= rd * r148;
      atomicAdd(vir, -(rd * r148));
      forcex = xx * r148;
      fxi += forcex;
      // f[j] -= forcex;
      atomicAdd(&(f[j]), -forcex);
      forcey = yy * r148;
      fyi += forcey;
      // f[j + 1] -= forcey;
      atomicAdd(&(f[j + 1]), -forcey);
      forcez = zz * r148;
      fzi += forcez;
      // f[j + 2] -= forcez;
      atomicAdd(&(f[j + 2]), -forcez);
    }
  }
  // f[i] += fxi;
  atomicAdd(&(f[i]), fxi);
  // f[i + 1] += fyi;
  atomicAdd(&(f[i + 1]), fyi);
  // f[i + 2] += fzi;
  atomicAdd(&(f[i + 2]), fzi);
}
