#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <omp.h>
#include "SML.h"


void BLOCK_INVERSE_ITERATION(BOX_BLOCK_II *Box_II){
   
   //Check input matrix
   if (Box_II->M->row_dim != Box_II->M->col_dim) {
      printf("Error in BLOCK_INVERSE_ITERATION\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box_II->M->row_dim, Box_II->M->col_dim);
      exit(1);
   }
   if (Box_II->M->row_dim <= 0 || Box_II->M->col_dim <= 0) {
      printf("Error in BLOCK_INVERSE_ITERATION\n");
      printf("row_dim(%d) or col_dim(%d) is illegal value\n", Box_II->M->row_dim, Box_II->M->col_dim);
      exit(1);
   }

   long i,j;
   int iter,conv_check;
   
   double *Error = GET_ARRAY_DOUBLE1(Box_II->eig_num);
   double *Un_Vector = GET_ARRAY_DOUBLE1(Box_II->M->row_dim);
   BOX_CG *Box_CG = malloc(sizeof(BOX_CG));
   
   FILE *file;
   
   for (iter = 0; iter <= Box_II->ii_max_step; iter++) {
      
      for (i = 0; i < Box_II->eig_num; i++){
         Error[i] = CONVERGE_CHECK(Box_II->M, Box_II->Eig_Vec[i], Box_II->Eig_Val[i], Box_II->p_threads);
      }
      
      conv_check = 0;
      for (i = 0; i < Box_II->eig_num; i++) {
         if (Error[i] < Box_II->ii_acc) {
            conv_check = conv_check + 1;
         }
      }
      
      //Check convergence
      if (conv_check == Box_II->eig_num) {
         
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/Block_inv_iter_step.txt","a+")) == NULL) {
            printf("Error in BLOCK_INVERSE_ITERATION\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         FREE_ARRAY_DOUBLE1(Un_Vector);
         FREE_ARRAY_DOUBLE1(Error);
         free(Box_CG);
         return;
      }
      
      for (i = 0; i < Box_II->eig_num; i++) {
         
         COPY_DOUBLE1(Box_II->Eig_Vec[i], Un_Vector, Box_II->M->row_dim, Box_II->p_threads);
         
         DIAG_ADD_CRS1(Box_II->M, Box_II->ii_diag_add - Box_II->Eig_Val[i], Box_II->p_threads);
         
         Box_CG->M         = Box_II->M;
         Box_CG->acc       = Box_II->cg_acc;
         Box_CG->Vec       = Box_II->Eig_Vec[i];
         Box_CG->Out_Vec   = Un_Vector;
         Box_CG->max_step  = Box_II->cg_max_step;
         Box_CG->p_threads = Box_II->p_threads;

         if (strcmp(Box_II->Type,"CG") == 0) {
            CONJUGATE_GRADIENT(Box_CG);
         }
         else if (strcmp(Box_II->Type,"MR") == 0) {
            MINIMUM_RESIDUAL(Box_CG);
         }
         else {
            printf("Error in BLOCL_INVERSE_ITERATION\n");
            printf("Invalid Type=%s\n", Box_II->Type);
            exit(1);
         }
         
         DIAG_ADD_CRS1(Box_II->M, -(Box_II->ii_diag_add - Box_II->Eig_Val[i]), Box_II->p_threads);
         
         NORMALIZE(Un_Vector, Box_II->M->row_dim, Box_II->p_threads);
         
#pragma omp parallel for num_threads (Box_II->p_threads)
         for (j = 0; j < Box_II->M->row_dim; j++) {
            Box_II->Eig_Vec[i][j] = Un_Vector[j];
         }
      }
      
      ORTHOGONALIZATION(Box_II->Eig_Vec, Box_II->eig_num, Box_II->M->row_dim, 0, Box_II->p_threads);
   }
   
   
   conv_check = 0;
   for (i = 0; i < Box_II->eig_num; i++) {
      if (Error[i] < Box_II->ii_acc) {
         conv_check = conv_check + 1;
      }
   }
   
   //Check convergence
   if (conv_check == Box_II->eig_num) {
      
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/Block_inv_iter_step.txt","a+")) == NULL) {
         printf("Error in BLOCK_INVERSE_ITERATION\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "%d\n", iter);
      fclose(file);
      FREE_ARRAY_DOUBLE1(Un_Vector);
      FREE_ARRAY_DOUBLE1(Error);
      free(Box_CG);
      return;
   }
   
   
   printf("Error in BLOCL_INVERSE_ITERATION\n");
   printf("Does not converge (max step is %d)\n", Box_II->ii_max_step);
   printf("continue...\n");
   
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file,"Error in BLOCK_INVERSE_ITERATION\n");
   fprintf(file,"Does not converge (max step is %d)\n",Box_II->ii_max_step);
   for (i = 0;i < Box_II->eig_num; i++) {
      double temp = CONVERGE_CHECK(Box_II->M, Box_II->Eig_Vec[i], Box_II->Eig_Val[i], Box_II->p_threads);
      printf("Error[%ld]=%e\n", i, temp);
      fprintf(file,"Error[%ld]=%e\n", i, temp);
   }
   fprintf(file,"continue...\n");
   fclose(file);
   
   FREE_ARRAY_DOUBLE1(Un_Vector);
   FREE_ARRAY_DOUBLE1(Error);
   free(Box_CG);
   
}

