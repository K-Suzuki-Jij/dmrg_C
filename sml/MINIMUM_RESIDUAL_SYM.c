#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <omp.h>
#include "SML.h"

void MINIMUM_RESIDUAL_SYM(BOX_CG *Box) {
   
   //Chech input matrix
   if (Box->M->row_dim != Box->M->col_dim) {
      printf("Error in MINIMUM_RESIDUAL_SYM\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box->M->row_dim, Box->M->col_dim);
      exit(1);
   }
   
   long i;
   int iter;
   int m_size = Box->M->row_dim;
   
   double *V0 = GET_ARRAY_DOUBLE1(m_size);
   double *V1 = GET_ARRAY_DOUBLE1(m_size);
   double *V2 = GET_ARRAY_DOUBLE1(m_size);
   double *W0 = GET_ARRAY_DOUBLE1(m_size);
   double *W1 = GET_ARRAY_DOUBLE1(m_size);
   double *W2 = GET_ARRAY_DOUBLE1(m_size);
   double *X1 = GET_ARRAY_DOUBLE1(m_size);
   double *Temp = GET_ARRAY_DOUBLE1(m_size);
   double **Temp_V = GET_ARRAY_DOUBLE2(Box->p_threads, m_size);
   double gamma0,gamma1,eta,delta,s0,s1,s2,c0,c1,c2,alpha0,alpha1,alpha2,alpha3,temp;
   FILE *file;


   MATRIX_VECTOR_PRODUCT_SYM(Box->M, Box->Out_Vec, V1, Temp_V, Box->p_threads);
   
#pragma omp parallel for num_threads (Box->p_threads)
   for (i = 0; i < m_size; i++) {
      V1[i] = Box->Vec[i] - V1[i];
   }
   
   
   gamma0 = L2_NORM(V1, m_size, Box->p_threads);
   eta = gamma0;
   s0 = 0;
   s1 = 0;
   c0 = 1;
   c1 = 1;
   
   for (iter = 0; iter < Box->max_step; iter++) {
      
      temp = 1.0/gamma0;
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V1[i] = V1[i]*temp;
      }
      
      MATRIX_VECTOR_PRODUCT_SYM(Box->M, V1, Temp, Temp_V, Box->p_threads);
      
      delta = INNER_PRODUCT(Temp, V1, m_size, Box->p_threads);
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V2[i] = Temp[i] - delta*V1[i] - gamma0*V0[i];
      }
      
      gamma1 = L2_NORM(V2, m_size, Box->p_threads);
      
      alpha0 = c1*delta - c0*s1*gamma0;
      alpha1 = sqrt(alpha0*alpha0 + gamma1*gamma1);
      alpha2 = s1*delta + c0*c1*gamma0;
      alpha3 = s0*gamma0;
      
      c2 = alpha0/alpha1;
      s2 = gamma1/alpha1;
      
      temp = 1.0/alpha1;
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         W2[i] = temp*(V1[i] - alpha3*W0[i] - alpha2*W1[i]);
         X1[i] = Box->Out_Vec[i] + c2*eta*W2[i];
      }

      eta = -s2*eta;
      fflush(stdout);
      printf("\rMINIMUM_RESIDUAL_error[%d]=%e", iter, fabs(eta));
      
      //Check convergence
      if (fabs(eta) < Box->acc) {
         
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/MINIMUM_RESIDUAL_Step.txt","a+")) == NULL) {
            printf("Error in MINIMUM_RESIDUAL_SYM\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         
         FREE_ARRAY_DOUBLE1(V0);
         FREE_ARRAY_DOUBLE1(V1);
         FREE_ARRAY_DOUBLE1(V2);
         FREE_ARRAY_DOUBLE1(W0);
         FREE_ARRAY_DOUBLE1(W1);
         FREE_ARRAY_DOUBLE1(W2);
         FREE_ARRAY_DOUBLE1(X1);
         FREE_ARRAY_DOUBLE1(Temp);
         FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);

         printf("\r                          ");
         printf("                          \r");
         return;
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < m_size; i++) {
         V0[i] = V1[i];
         V1[i] = V2[i];
         
         W0[i] = W1[i];
         W1[i] = W2[i];
         
         Box->Out_Vec[i] = X1[i];
      }
      
      s0 = s1;
      s1 = s2;
      
      c0 = c1;
      c1 = c2;
      
      gamma0 = gamma1;
      
   }
   
   printf("Error in MINIMUM_RESIDUAL_SYM\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   printf("Continue...\n");
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in MINIMUM_RESIDUAL_SYM\n");
   fprintf(file, "Does not converge (max step is %d)\n", Box->max_step);
   fprintf(file, "Continue...\n");
   fclose(file);
   
   FREE_ARRAY_DOUBLE1(V0);
   FREE_ARRAY_DOUBLE1(V1);
   FREE_ARRAY_DOUBLE1(V2);
   FREE_ARRAY_DOUBLE1(W0);
   FREE_ARRAY_DOUBLE1(W1);
   FREE_ARRAY_DOUBLE1(W2);
   FREE_ARRAY_DOUBLE1(X1);
   FREE_ARRAY_DOUBLE1(Temp);
   FREE_ARRAY_DOUBLE2(Temp_V, Box->p_threads);
   
}
