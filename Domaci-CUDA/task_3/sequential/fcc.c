
/*
 *   Generate fcc lattice for atoms inside the box
 */
  void
  fcc(double x[], int npart, int mm, double a){
    int ijk=0;
    int i,j,k,lg;
  //ovde moze lagana sekcija, ali da se izracuna ijk u donjem delu, ugl podaci u nizu ne zavise od prethodnih nego samo od svoje pozicije
  //popunjava se prva polovina niza x(:,:,:,0:1)
    for (lg=0; lg<2; lg++)
      for (i=0; i<mm; i++)
        for (j=0; j<mm; j++)
          for (k=0; k<mm; k++) {
            x[ijk]   = i*a+lg*a*0.5;
            x[ijk+1] = j*a+lg*a*0.5;
            x[ijk+2] = k*a;
            ijk += 3;
          }
  //popunjava se prva polovina niza x(:,:,:,2:3)
    for (lg=1; lg<3; lg++)
      for (i=0; i<mm; i++)
        for (j=0; j<mm; j++)
          for (k=0; k<mm; k++) {
            x[ijk]   = i*a+(2-lg)*a*0.5;
            x[ijk+1] = j*a+(lg-1)*a*0.5;
            x[ijk+2] = k*a+a*0.5;
            ijk += 3;
          }

  }
