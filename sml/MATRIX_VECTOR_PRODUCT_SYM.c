#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"
#include <omp.h>

void MATRIX_VECTOR_PRODUCT_SYM(CRS1 *M, double *V1, double *Out, double **Temp, int p_threads) {
   
   //Note that all the diagonal elements must be stored even if they are zero.
   
   int thread_num;
   long i,j;
   double temp1,temp2;
   
   int start_row = 0;
   while (M->Row[start_row+1] - 1 < 0) {
      start_row++;
   }
   
   if (start_row != 0) {
      printf("Error in MATRIX_VECTOR_PRODUCT_SYM\n");
      printf("start_row=%d\n",start_row);
      exit(1);
   }

#pragma omp parallel for num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      Out[i] = 0;
   }
   
#pragma omp parallel for private (j,thread_num,temp1,temp2) schedule(guided) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      thread_num = omp_get_thread_num();
      temp1      = M->Val[M->Row[i+1] - 1]*V1[i];
      temp2      = V1[i];
      for (j = M->Row[i]; j < M->Row[i+1] - 1; j++) {
         temp1 = temp1 + M->Val[j]*V1[M->Col[j]];
         Temp[thread_num][M->Col[j]] = Temp[thread_num][M->Col[j]] +  M->Val[j]*temp2;
      }
      Temp[thread_num][i] = Temp[thread_num][i] + temp1;
   }
   
#pragma omp parallel for private (j) num_threads (p_threads)
   for (i = 0; i < M->row_dim; i++) {
      for (j = 0; j < p_threads; j++) {
         Out[i] = Out[i] + Temp[j][i];
         Temp[j][i] = 0;
      }
   }
    
}
