// create coarse graining time series from the original time series.
// the return pointer need to be deleted after being used.
#include "helper.h"
double* coarse_grain(double* data, int len, int s)
{
  int cg_len = len/s;
  double* cg_data = new double[cg_len]();
  
  for (int i = 0; i < cg_len; i++) {
    for (int k = 0; k < s; k++)
      cg_data[i] += data[i*s+k];
    cg_data[i] /= s;
  }
  
  return cg_data;
}

int cmpfunc (const void * a, const void * b)
{
  return ((struct Item *)a)->value > ((struct Item *)b)->value ? 1 : -1;
}
