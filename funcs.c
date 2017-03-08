#include "funcs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Func1(unsigned char a) {
  double x = a;
  x *= 0.001;
  double y = sqrt(x) + 0.12;
  double z = pow(x, 1.3);
  double w = z / y + 0.4;
  return pow(w, 2) / 25;
}

double Func2(unsigned char a, unsigned char b) {
  int aa = a;
  int bb = b;
  int mid = (aa + bb) / 2;
  int step = 1;
  if(bb < aa)
    step = -1;
  int k;
  double sum = 0;
  for(k = aa; k != mid; k += step)
    sum += Func1(k) / 300;
  return sum;
}

// This function expects the input array a to consist of 4 numbers.
double Func3(unsigned char* a, double p) {
  double* v = malloc(4*sizeof(double));
  double* w = malloc(4*sizeof(double));
  int k;
  for(k = 0; k < 4; k++) {
    v[k] = a[k];
    v[k] *= 0.001;
  }
  w[0] = v[0] - 0.6*v[1] + 0.4*v[2];
  w[1] = v[1] - 0.5*v[2] + 0.9*v[3];
  w[2] = v[0] - 0.4*v[2] + 0.1*v[3];
  w[3] = v[0] - 0.2*v[1] + 0.3*v[2];
  double sum = 0;
  for(k = 0; k < 4; k++)
    sum += pow(w[k], p) / 3;
  free(v);
  free(w);
  return sum;
}

double ComputeNumber(unsigned char* buf, int nBytes, double p) {
  int i;
  double sum = 0;
  for(i = 0; i < nBytes; i++)
    sum += Func1(buf[i]);
  for(i = 0; i < nBytes-1; i+=2)
    sum += Func2(buf[i], buf[i+1]);
  for(i = 0; i < nBytes-3; i+=4)
    sum += Func3(&buf[i], p);
  return sum;
}
