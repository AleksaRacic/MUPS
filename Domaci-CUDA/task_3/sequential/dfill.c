
/*
 *  function dfill : intialises double precision array to scalar value
 */
  void
  dfill(int n, double val, double a[], int ia){
    int i;
    //ovo je jednopstavna paralelizacija
    for (i=0; i<(n-1)*ia+1; i+=ia)
      a[i] = val;
  }
