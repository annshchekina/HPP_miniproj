#include "funcs.h"
#include "timings.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <immintrin.h>

#define characteristic_coeff 3
#define VECT_LEN 4

double sqrt_Func1, pow_Func1;
const double inv_Func2 = 1 / 300, inv_Func3 = 1 / 3;

__m256d a_Func1, b_Func1, c_Func1, d_Func1, e_Func1; 

double Func1(unsigned char a, double * time_sec) 
{ 
  double start = get_wall_seconds();
  double y = sqrtf(a) * sqrt_Func1 + 0.12;
  double z = powf(a, 1.3) * pow_Func1;
  double w = z / y + 0.4;
  double x = w * w * 0.04;
  *time_sec += get_wall_seconds() - start;
  return x;
}

__m256d Func1_avx(const __m256d v1, const __m256d v2, 
  double * time_sec) 
{  
  double start = get_wall_seconds();
  __m256d x = _mm256_mul_pd(v2, c_Func1);
  __m256d y = _mm256_add_pd(_mm256_mul_pd(v1, a_Func1), b_Func1);
  __m256d z = _mm256_div_pd(x, y);
  __m256d w = _mm256_add_pd(z, d_Func1);
  __m256d r = _mm256_mul_pd(_mm256_mul_pd(w, w), e_Func1);
  *time_sec += get_wall_seconds() - start;
  return r;
}

double Func2(unsigned char a, unsigned char b, 
  double * time_sec, double * time_sec_Func1) 
{
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
double Func3(unsigned char* a, double p, double * time_sec) 
{
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
	double * Func2_time, double * Func3_time) 
{
  assert(p == characteristic_coeff);
  
  sqrt_Func1 = sqrt(0.001);
  pow_Func1 = pow(0.001, 1.3);
  
  a_Func1 = _mm256_set1_pd(sqrt_Func1);
  b_Func1 = _mm256_set1_pd(0.12);
  c_Func1 = _mm256_set1_pd(pow_Func1);
  d_Func1 = _mm256_set1_pd(0.4);
  e_Func1 = _mm256_set1_pd(0.04);
  
  int i, j;
  double * v1_Func1 = (double *)malloc(nBytes * sizeof(double)), 
  	* v2_Func1 = (double *)malloc(nBytes * sizeof(double));
  
  double start = get_wall_seconds();
  for(i = 0; i < nBytes; i++)
  {
	  
	  v1_Func1[i] = sqrtf(buf[i]);
	  v2_Func1[i] = powf(buf[i], 1.3);
  }
  *Func1_time += get_wall_seconds() - start;
  
  double sum = 0;
 
  int n = nBytes - nBytes % VECT_LEN;
  for(i = 0; i < n; i += VECT_LEN)
  { 
	__m256d v1 = _mm256_loadu_pd(v1_Func1 + i); 
	__m256d v2 = _mm256_loadu_pd(v2_Func1 + i); 
	__m256d v = Func1_avx(v1, v2, Func1_time);
	
	double res_vec[VECT_LEN];
	_mm256_storeu_pd(res_vec, v);
	for(j = 0; j < VECT_LEN; j++)
		sum += res_vec[j];
  }
  for(i = n; i < nBytes; i++)
  {
	  sum += Func1(buf[i], Func1_time);
  }
  
  /*for(i = 0; i < nBytes; i++)
  {
  	sum += Func1(buf[i], Func1_time);
  }*/
  
  for(i = 0; i < nBytes - 1; i += 2)
    sum += Func2(buf[i], buf[i+1], Func2_time, Func1_Func2_time);
  
  for(i = 0; i < nBytes - 3; i += 4)
    sum += Func3(&buf[i], p, Func3_time);
    
  free(v1_Func1);
  free(v2_Func1);
  
  return sum;
}
