#include <stdio.h>
#include <stdlib.h>
#include <math.h>

  void srand48(long);
  double drand48(void);
/*
 *  Sample Maxwell distribution at temperature tref
 */
  void
  mxwell(double vh[], int n3, double h, double tref){ //n3 je npart * 3
    int i;
    int npart=n3/3;
    double r, tscale, v1, v2, s, ekin=0.0, sp=0.0, sc;
    
    srand48(4711);
    tscale=16.0/((double)npart-1.0);


    //mozed paralelizacija
    for (i=0; i<n3; i+=2) {
      s=2.0;
      while (s>=1.0) {
        v1=2.0*drand48()-1.0;
        v2=2.0*drand48()-1.0;
        s=v1*v1+v2*v2;
      }
      r=sqrt(-2.0*log(s)/s); // r moze biti lokalno za svaku nit
      vh[i]=v1*r;
      vh[i+1]=v2*r;
    }
    

    //IZRACUNAVANJE SP_X, SP_Y, SP_Z odnosno sp[3] se moze staviti u onu gore for petlju s[i % 3] += vh
    //Onda treba proci kroz ceo niz i oduzeti odgovarajuci sp,  aredukcija je po ekin

    //ovde se radi normalizacija(da mean bude 0 po x osi)
    for (i=0; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=0; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }
    //isto po y osi
    sp=0.0;
    for (i=1; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=1; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }

    //isto po z osi
    sp=0.0;
    for (i=2; i<n3; i+=3) sp+=vh[i];
    sp/=(double)npart;
    for(i=2; i<n3; i+=3) {
      vh[i]-=sp;
      ekin+=vh[i]*vh[i];
    }

    sc=h*sqrt(tref/(tscale*ekin));
    for (i=0; i<n3; i++) vh[i]*=sc;
  }
