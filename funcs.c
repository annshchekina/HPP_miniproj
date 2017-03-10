#include "funcs.h"
#include "timings.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define characteristic_coeff 3

double sqrt_Func1, pow_Func1;
const double inv_Func2 = 1 / 300, inv_Func3 = 1 / 3;

double Func1(unsigned char a, double * time_sec) { // 1B is a char
  double start = get_wall_seconds();
  double y = sqrtf(a) * sqrt_Func1 + 0.12;
  double z = powf(a, 1.3) * pow_Func1;
  double w = z / y + 0.4;
  double x = w * w * 0.04;
  *time_sec += get_wall_seconds() - start;
  return x;
}

double Func2(unsigned char a, unsigned char b, double * time_sec, double * time_sec_Func1) {
  double start = get_wall_seconds();
  int aa = a;
  int bb = b;
  int mid = (aa + bb) / 2;
  int step = 1;
  if(bb < aa)
    step = -1;
  int k;
  double sum = 0;
  for(k = aa; k != mid; k += step)
    sum += Func1(k, time_sec_Func1);
  sum *= inv_Func2;
  *time_sec += get_wall_seconds() - start;
  return sum;
}

// This function expects the input array a to consist of 4 numbers.
double Func3(unsigned char* a, double p, double * time_sec) {
  double start = get_wall_seconds();
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
    sum += w[k] * w[k] * w[k]; // power characteristic_coeff
  free(v);
  free(w);
  sum *= inv_Func3;
  *time_sec += get_wall_seconds() - start;
  return sum;
}

double ComputeNumber(unsigned char* buf, int nBytes, double p, 
	double * Func1_time, double * Func1_Func2_time, 
	double * Func2_time, double * Func3_time) {
  assert(p == characteristic_coeff);
  sqrt_Func1 = sqrt(0.001);
  pow_Func1 = pow(0.001, 1.3);
  int i;
  double sum = 0;
  for(i = 0; i < nBytes; i++)
    sum += Func1(buf[i], Func1_time);
  for(i = 0; i < nBytes-1; i+=2)
    sum += Func2(buf[i], buf[i+1], Func2_time, Func1_Func2_time);
  for(i = 0; i < nBytes-3; i+=4)
    sum += Func3(&buf[i], p, Func3_time);
  return sum;
}
