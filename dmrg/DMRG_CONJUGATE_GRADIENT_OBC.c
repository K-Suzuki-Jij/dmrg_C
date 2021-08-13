#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"
#include "dmrg.h"

void DMRG_CONJUGATE_GRADIENT_OBC(DMRG_BOX_CG *Box) {
   
   long i;
   int iter;
   double alpha,beta,inner_pro_rrr,temp;
   int dim_LLLRRRRL = Box->Dmrg_Basis->dim_LLLRRRRL;
   double diag_val = Box->diag_val;
   double *RRR = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   double *PPP = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   double *YYY = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   FILE *file;

   //Initial vector Box->Out_Vec must be a nonzero vector
   //Normalize initial vector
   NORMALIZE(Box->Out_Vec, dim_LLLRRRRL, Box->p_threads);
   
   DMRG_MATRIX_VECTOR_PRODUCT_OBC(Box->M_LLLR, Box->M_LRRL, Box->M_LRRL_Sign, Box->M_RRRL, Box->Out_Vec, RRR, Box->Ele_RR, Box->p_threads, Box->Dmrg_Basis);
   
#pragma omp parallel for num_threads (Box->p_threads)
   for (i = 0; i < dim_LLLRRRRL; i++) {
      RRR[i] = Box->Vec[i] - RRR[i] - Box->Out_Vec[i]*diag_val;
      PPP[i] = RRR[i];
   }
   
   
   for (iter = 1; iter < Box->max_step; iter++) {
      
      DMRG_MATRIX_VECTOR_PRODUCT_OBC(Box->M_LLLR, Box->M_LRRL, Box->M_LRRL_Sign, Box->M_RRRL, PPP, YYY, Box->Ele_RR, Box->p_threads, Box->Dmrg_Basis);
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < dim_LLLRRRRL; i++) {
         YYY[i] = YYY[i] + PPP[i]*diag_val;
      }
      
      alpha = INNER_PRODUCT(RRR, RRR, dim_LLLRRRRL, Box->p_threads);
      inner_pro_rrr = alpha;
      alpha = alpha/INNER_PRODUCT(PPP, YYY, dim_LLLRRRRL, Box->p_threads);
    
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < dim_LLLRRRRL; i++) {
         Box->Out_Vec[i] = Box->Out_Vec[i] + alpha*PPP[i];
      }
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < dim_LLLRRRRL; i++) {
         RRR[i] = RRR[i] - alpha*YYY[i];
      }
      
      temp = INNER_PRODUCT(RRR, RRR, dim_LLLRRRRL, Box->p_threads);
      
      fflush(stdout);
      printf("\rCONJUGATE_GRADIENT_error[%d]=%e",iter, temp);

      //Check convergence
      if (temp < Box->acc) {
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/CONJUGATE_GRADIENT_Step.txt","a+")) == NULL) {
            printf("Error in DMRG_CONJUGATE_GRADIENT_OBC\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         
         FREE_ARRAY_DOUBLE1(RRR);
         FREE_ARRAY_DOUBLE1(PPP);
         FREE_ARRAY_DOUBLE1(YYY);
         printf("\r                              ");
         printf("                              \r");
      
         return;
      }

      beta = temp/inner_pro_rrr;
      
#pragma omp parallel for num_threads (Box->p_threads)
      for (i = 0; i < dim_LLLRRRRL; i++) {
         PPP[i] = RRR[i] + beta*PPP[i];
      }
   }
   
   printf("Error in DMRG_CONJUGATE_GRADIENT_OBC\n");
   printf("Does not converge (max step is %d)\n", Box->max_step);
   printf("Continue...\n");
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in DMRG_CONJUGATE_GRADIENT_OBC\n");
   fprintf(file, "Does not converge (max step is %d)\n" ,Box->max_step);
   fprintf(file, "Continue...\n");
   fclose(file);
   if ((file = fopen("./SML_out/CONJUGATE_GRADIENT_Step.txt","a+")) == NULL) {
      printf("Error in DMRG_CONJUGATE_GRADIENT_OBC\n");
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "%d\n", Box->max_step);
   fclose(file);
   
   FREE_ARRAY_DOUBLE1(RRR);
   FREE_ARRAY_DOUBLE1(PPP);
   FREE_ARRAY_DOUBLE1(YYY);
   
}
