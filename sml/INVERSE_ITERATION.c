#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"

void INVERSE_ITERATION(BOX_II *Box_II) {
   
   //Check input matrix
   if (Box_II->M->row_dim != Box_II->M->col_dim) {
      printf("Error in INVERSE_ITERATION\n");
      printf("The input matrix is not a square matirx(row=%d,col=%d)\n", Box_II->M->row_dim, Box_II->M->col_dim);
      exit(1);
   }
   
   double *Un_Vector = GET_ARRAY_DOUBLE1(Box_II->M->row_dim);
   double error      = 1000;
   BOX_CG *Box_CG    = malloc(sizeof(BOX_CG));
   int iter;
   FILE *file;
   
   DIAG_ADD_CRS1(Box_II->M, Box_II->ii_diag_add - Box_II->Eig_Val[0], Box_II->p_threads);
   
   COPY_DOUBLE1(Box_II->Eig_Vec, Un_Vector, Box_II->M->row_dim, Box_II->p_threads);
   
   for (iter = 0; iter <= Box_II->ii_max_step; iter++) {
      
      error = CONVERGE_CHECK(Box_II->M, Box_II->Eig_Vec, Box_II->ii_diag_add, Box_II->p_threads);
    
      //Check convergence
      if (error < Box_II->ii_acc) {
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/Inv_iter_step.txt","a+")) == NULL) {
            printf("Error in INVERSE_ITERATION\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         
         Box_II->error = error;
         FREE_ARRAY_DOUBLE1(Un_Vector);
         
         DIAG_ADD_CRS1(Box_II->M, -(Box_II->ii_diag_add - Box_II->Eig_Val[0]), Box_II->p_threads);
         
         free(Box_CG);
         
         return;
      }
      
      Box_CG->M         = Box_II->M;
      Box_CG->acc       = Box_II->cg_acc;
      Box_CG->Vec       = Box_II->Eig_Vec;
      Box_CG->Out_Vec   = Un_Vector;
      Box_CG->max_step  = Box_II->cg_max_step;
      Box_CG->p_threads = Box_II->p_threads;

      if (strcmp(Box_II->Type,"CG") == 0) {
         CONJUGATE_GRADIENT(Box_CG);
      }
      else if(strcmp(Box_II->Type,"MR") == 0){
         MINIMUM_RESIDUAL(Box_CG);
      }
      else {
         printf("Error in INVERSE_ITERATION\n");
         printf("Invalid Box_II->Type=%s\n",Box_II->Type);
         exit(1);
      }
      
      NORMALIZE(Un_Vector, Box_II->M->row_dim, Box_II->p_threads);
      
      COPY_DOUBLE1(Un_Vector, Box_II->Eig_Vec, Box_II->M->row_dim, Box_II->p_threads);
      
   }
   
   
   //Check convergence
   if (error < Box_II->ii_acc) {
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/Inv_iter_step.txt","a+")) == NULL) {
         printf("Error in INVERSE_ITERATION\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "%d\n", iter);
      fclose(file);
      Box_II->error = error;
      FREE_ARRAY_DOUBLE1(Un_Vector);
      
      DIAG_ADD_CRS1(Box_II->M, -(Box_II->ii_diag_add - Box_II->Eig_Val[0]), Box_II->p_threads);
      
      free(Box_CG);
      
      return;
   }
   
   printf("Error in INVERSE_ITERATION\n");
   printf("Does not converge(%e > %e) (max step is %d)\n", error, Box_II->ii_acc ,Box_II->ii_max_step);
   printf("continue...\n");
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in INVERSE_ITERATION\n");
   fprintf(file, "Does not converge(%e > %e) (max step is %d)\n", error, Box_II->ii_acc, Box_II->ii_max_step);
   fprintf(file, "continue...\n");
   fclose(file);
   
   DIAG_ADD_CRS1(Box_II->M, -(Box_II->ii_diag_add - Box_II->Eig_Val[0]), Box_II->p_threads);
   FREE_ARRAY_DOUBLE1(Un_Vector);
   free(Box_CG);
   
}

