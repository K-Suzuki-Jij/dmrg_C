#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "SML.h"

double CONVERGE_CHECK(CRS1 *M, double *Eig_Vec, double eig_val, int p_threads) {
   
   long i,j;
   double temp,norm;
   norm = 0;
   
#pragma omp parallel for private (j,temp) reduction (+:norm) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      temp = 0;
      for (j = M->Row[i]; j < M->Row[i+1]; j++) {
         temp = temp + M->Val[j]*Eig_Vec[M->Col[j]];
      }
      norm = norm + fabs(temp - eig_val*Eig_Vec[i]);
   }
   
   return norm;
}
