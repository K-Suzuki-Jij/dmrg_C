#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SML.h"

void ORTHOGONALIZATION(double **Vector, long v_num, long dim, long start, int p_threads) {
   
   if (v_num == 1) {
      NORMALIZE(Vector[0], dim, p_threads);
      return;
   }
   
   long i,j,k;
   double *Temp = GET_ARRAY_DOUBLE1(dim);
   
   for (i = start; i < v_num; i++) {
      
#pragma omp parallel for num_threads (p_threads)
      for (j = 0; j < dim; j++) {
         Temp[j] = Vector[i][j];
      }
      
      for (j = 0; j < i; j++) {
         double inn_pro = INNER_PRODUCT(Vector[j], Temp, dim, p_threads);
         
#pragma omp parallel for num_threads (p_threads)
         for (k = 0; k < dim; k++) {
            Temp[k] = Temp[k] - inn_pro*Vector[j][k];
         }
      }
      
      NORMALIZE(Temp, dim, p_threads);
      
#pragma omp parallel for num_threads (p_threads)
      for (j = 0; j < dim; j++) {
         Vector[i][j] = Temp[j];
      }
   }
   
   FREE_ARRAY_DOUBLE1(Temp);
}


