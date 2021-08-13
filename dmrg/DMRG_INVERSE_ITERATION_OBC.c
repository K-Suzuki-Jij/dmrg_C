#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/stat.h>
#include "SML.h"
#include "dmrg.h"

void DMRG_INVERSE_ITERATION_OBC(DMRG_BOX_II *Box_II) {
   
   int dim_LLLRRRRL  = Box_II->Dmrg_Basis->dim_LLLRRRRL;
   double *Un_Vector = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   double error      = 1000;
   double eig_val    = Box_II->eig_val[0];
   double diag_val   = Box_II->ii_diag_add - eig_val;
   int iter;
   FILE *file;
   DMRG_BOX_CG *Box_CG = malloc(sizeof(DMRG_BOX_CG));

   for (iter = 0; iter <= Box_II->ii_max_step; iter++) {
      
      error = DMRG_CONVERGE_CHECK_OBC(eig_val, Box_II->M_LLLR, Box_II->M_LRRL, Box_II->M_LRRL_Sign, Box_II->M_RRRL, Box_II->Eig_Vec, Un_Vector, Box_II->Ele_RR, Box_II->p_threads, Box_II->Dmrg_Basis);
    
      //Check convergence
      if (error < Box_II->ii_acc) {
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/Inv_iter_step.txt","a+")) == NULL) {
            printf("Error in DMRG_INVERSE_ITERATION_OBC\n");
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file, "%d\n", iter);
         fclose(file);
         
         Box_II->error = error;
         
         FREE_ARRAY_DOUBLE1(Un_Vector);
         
         free(Box_CG);
         
         return;
      }
      
      COPY_DOUBLE1(Box_II->Eig_Vec, Un_Vector, dim_LLLRRRRL, Box_II->p_threads);
      
      Box_CG->acc         = Box_II->cg_acc;
      Box_CG->diag_val    = diag_val;
      Box_CG->Vec         = Box_II->Eig_Vec;
      Box_CG->Out_Vec     = Un_Vector;
      Box_CG->max_step    = Box_II->cg_max_step;
      Box_CG->p_threads   = Box_II->p_threads;
      Box_CG->M_LLLR      = Box_II->M_LLLR;
      Box_CG->M_LRRL      = Box_II->M_LRRL;
      Box_CG->M_LRRL_Sign = Box_II->M_LRRL_Sign;
      Box_CG->M_RRRL      = Box_II->M_RRRL;
      Box_CG->Ele_RR      = Box_II->Ele_RR;
      Box_CG->Dmrg_Basis  = Box_II->Dmrg_Basis;

      DMRG_CONJUGATE_GRADIENT_OBC(Box_CG);
      
      NORMALIZE(Un_Vector, dim_LLLRRRRL, Box_II->p_threads);
      
      COPY_DOUBLE1(Un_Vector, Box_II->Eig_Vec, dim_LLLRRRRL, Box_II->p_threads);
      
   }
   
   //Check convergence
   if (error < Box_II->ii_acc) {
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/Inv_iter_step.txt","a+")) == NULL) {
         printf("Error in DMRG_INVERSE_ITERATION_OBC\n");
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file, "%d\n", iter);
      fclose(file);
      
      Box_II->error = error;
      
      FREE_ARRAY_DOUBLE1(Un_Vector);
            
      free(Box_CG);
      
      return;
   }
   
   printf("Error in DMRG_INVERSE_ITERATION_OBC\n");
   printf("Does not converge(%e > %e) (max step is %d)\n", error, Box_II->ii_acc ,Box_II->ii_max_step);
   printf("continue...\n");
   mkdir("SML_out",0777);
   if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
      printf("Can't open file.\n");
      exit(1);
   }
   fprintf(file, "Error in DMRG_INVERSE_ITERATION_OBC\n");
   fprintf(file, "Does not converge(%e > %e) (max step is %d)\n", error, Box_II->ii_acc, Box_II->ii_max_step);
   fprintf(file, "continue...\n");
   fclose(file);
   
   Box_II->error = error;
   FREE_ARRAY_DOUBLE1(Un_Vector);
   free(Box_CG);
   
}

