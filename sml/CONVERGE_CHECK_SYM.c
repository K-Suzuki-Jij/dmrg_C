#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "SML.h"

double CONVERGE_CHECK_SYM(CRS1 *M, double *Eig_Vec, double eig_val, int p_threads) {
   
   long i,j;
   double temp1,temp2,norm;
   
   int start_row = 0;
   while (M->Row[start_row+1] - 1 < 0) {
      start_row++;
   }
   
   if (start_row != 0) {
      printf("Error in CONVERGE_CHECK_SYM\n");
      exit(1);
   }
   
   double *Temp_V = GET_ARRAY_DOUBLE1(M->row_dim);
      
   for (i = 0; i < M->row_dim; i++) {
      temp1 = M->Val[M->Row[i + 1] - 1]*Eig_Vec[i];
      temp2 = Eig_Vec[i];
      for (j = M->Row[i]; j < M->Row[i+1] - 1; j++) {
         temp1 = temp1 + M->Val[j]*Eig_Vec[M->Col[j]];
         Temp_V[M->Col[j]] = Temp_V[M->Col[j]] +  M->Val[j]*temp2;
      }
      Temp_V[i] = Temp_V[i] + temp1;
   }
   
   norm = 0;
#pragma omp parallel for reduction (+:norm) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      norm = norm + fabs(eig_val*Eig_Vec[i] - Temp_V[i]);
   }
   
   FREE_ARRAY_DOUBLE1(Temp_V);

   
   return norm;
}
