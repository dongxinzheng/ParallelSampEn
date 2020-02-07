#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

// declarations
double do_basic_sample_entropy(double* data, int len, int m, double r, double sd);
double* coarse_grain(double* data, int len, int s);

/*
############################################################################################
# Function to calculate sample entropy using Costa's MSE c code (2004)
#
# ARGUMENTS
#  data: A numeric vector representing a time series.
#  m: A integer value, length of template vector.
#  r: A double value, coefficient used to calculate similarity thresholds by r*sd(data)
#  s: A integer value, scale factor for coarse graining of the time series.
# VALUE
#  the sample entroy value of specified parameters: m, r, s.
############################################################################################
*/

// [[Rcpp::export]]
NumericVector basic_sample_entropy(NumericVector data, IntegerVector m, NumericVector r, IntegerVector s) {
  double sd_value = sd(noNA(data));
  double* cg_data = coarse_grain(data.begin(),data.length(),s[0]);
  double se = do_basic_sample_entropy(cg_data,data.length()/s[0],m[0],r[0],sd_value);
  delete[] cg_data; 
  return NumericVector::create(se);
}

// internal function
double do_basic_sample_entropy(double* data, int len, int m, double r, double sd)
{
  int i,j;
  vector<long long> cont(m+2,0);
  double r_new = r*sd;
  int vec_count = len - m;
  
  for (i = 0; i < vec_count; ++i) {
    for (j = i+1; j < vec_count; ++j) {
      int k = 0;
      while (k <= m && fabs(data[i+k] - data[j+k]) <= r_new)
        cont[++k]++;
    }
  }
  
  double se;
  if (cont[m+1] == 0 || cont[m] == 0)
    se = -log((double)1/((vec_count)*(vec_count-1)));
  else
    se = -log((double)cont[m+1]/cont[m]);
  
  return se;
}


