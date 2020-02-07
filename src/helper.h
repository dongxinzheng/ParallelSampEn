#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

struct Item
{
  double value;
  int pos;
  bool operator < (const Item &m)const
  {
    return value < m.value;
  }
};
int cmpfunc (const void * a, const void * b);

double* coarse_grain(double* data, int len, int s);

#endif // HELPER_H_INCLUDED