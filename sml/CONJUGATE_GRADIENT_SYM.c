#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"

void CONJUGATE_GRADIENT_SYM(BOX_CG *Box){
   
   //Chech input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in CONJUGATE_GRADIENT_SYM\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n",Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   if (Box->M->row_dim <= 0 || Box->M->col_dim <= 0) {
      printf("Error in CONJUGATE_GRADIENT_SYM\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n",Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   
   long i;
   int iter;
   int m_size = Box->M->row_dim;
   double alpha,beta,inner_pro_rrr,temp;
   double *RRR = GET_ARRAY_DOUBLE1(m_size);
   double *PPP = GET_ARRAY_DOUBLE1(m_size);
   double *YYY = GET_ARRAY_DOUBLE1(m_size);
   double **Temp_V = GET_ARRAY_DOUBLE2(Box->p_threads, m_size);
   
   FILE *file;
   
   //Initial vector Box->Out_Vec must be a nonzero vector
   //Normalize initial vector
   NORMALIZE(Box->Out_Vec, m_size, Box->p_threads);
   
   MATRIX_VECTOR_PRODUCT_SYM(Box->M, Box->Out_Vec, RRR, Temp_V, Box->p_threads);
   
#pragma omp parallel for num_threads (Box->p_threads)
   for (i = 0; i < m_size; i++) {
      RRR[i] = Box->Vec[i] - RRR[i];
      PPP[i] = RRR[i];
   }
   
   for (iter = 1; iter < Box->max_step; iter++) {
      
      MATRIX_VECTOR_PRODUCT_SYM(Box->M, PPP, YYY, Temp_V, Box->p_threads);
      
      alpha = INNER_PRODUCT(RRR, RRR, m_size, Box->p_threads);
      inner_pro_rrr = alpha;
      alpha = alpha/INNER_PRODUCT(PPP, YYY, m_size, Box->p_threads);
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         Box->Out_Vec[i] = Box->Out_Vec[i] + alpha*PPP[i];
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         RRR[i] = RRR[i] - alpha*YYY[i];
      }
      
      temp = INNER_PRODUCT(RRR, RRR, m_size, Box->p_threads);
      
      fflush(stdout);
      printf("\rCONJUGATE_GRADIENT_error[%d]=%e",iter, temp);
      
      //Check convergence
      if (temp < Box->acc) {
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/CONJUGATE_GRADIENT_SYM_Step.txt","a+")) == NULL) {
            printf("Error in CONJUGATE_GRADIENT_SYM\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         
         FREE_ARRAY_DOUBLE1(RRR);
         FREE_ARRAY_DOUBLE1(PPP);
         FREE_ARRAY_DOUBLE1(YYY);
         FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);
         
         printf("\r                              ");
         printf("                              \r");
         
         return;
      }
      
      beta = temp/inner_pro_rrr;
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         PPP[i] = RRR[i] + beta*PPP[i];
      }
   }
   
   printf("Error in CONJUGATE_GRADIENT_SYM\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   printf("Continue...\n");
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in CONJUGATE_GRADIENT_SYM\n");
   fprintf(file, "Does not converge (max step is %d)\n" ,Box->max_step);
   fprintf(file, "Continue...\n");
   fclose(file);
   if ((file = fopen("./SML_out/CONJUGATE_GRADIENT_SYM_Step.txt","a+")) == NULL) {
      printf("Error in CONJUGATE_GRADIENT_SYM\n");
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "%d\n", Box->max_step);
   fclose(file);
   
   FREE_ARRAY_DOUBLE1(RRR);
   FREE_ARRAY_DOUBLE1(PPP);
   FREE_ARRAY_DOUBLE1(YYY);
   FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);
   
}
