#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "exact.h"

void EXACT_DIAGONALIZE_HAMILTONIAN(EXACT_HAM_INFO *Ham_Info, EXACT_PARAMETER *Param, EXACT_TIME *Time, int p_threads) {
   
   Time->diag = omp_get_wtime();
   
   Ham_Info->diag_num = Param->diag_num;

   BOX_LAN *Box_Lan = malloc(sizeof(BOX_LAN));
   Box_Lan->M              = Ham_Info->Ham;
   Box_Lan->Eig_Vec        = Ham_Info->Vector[0];
   Box_Lan->eig_val        = &Ham_Info->Value[0];
   Box_Lan->acc            = Param->diag_acc;
   Box_Lan->min_step       = Param->diag_min_step;
   Box_Lan->max_step       = Param->diag_max_step;
   Box_Lan->max_block_step = Param->diag_max_step;
   Box_Lan->p_threads      = p_threads;
   strcpy(Box_Lan->Type, "Normal");
   strcpy(Box_Lan->Guess,"No");
   
   if (strcmp(Param->Diag_Method, "Lanczos") == 0) {
      LANCZOS(Box_Lan);
   }
   else if (strcmp(Param->Diag_Method, "Lanczos_Slow") == 0) {
      LANCZOS_SLOW(Box_Lan);
   }
   else {
      printf("Error in DMRG_DIAGONALIZE_HAMILTONIAN\n");
      exit(1);
   }
   
   free(Box_Lan);
   
   Time->diag = omp_get_wtime() - Time->diag;
   
   //Inverse Iteration
   Time->inv_iter = omp_get_wtime();
   
   BOX_II *Box_II = malloc(sizeof(BOX_II));
   Box_II->M           = Ham_Info->Ham;
   Box_II->ii_acc      = Param->inv_iter_acc;
   Box_II->ii_diag_add = Param->inv_iter_diag_add;
   Box_II->ii_max_step = Param->inv_iter_max_step;
   Box_II->cg_acc      = Param->cg_acc;
   Box_II->cg_max_step = Param->cg_max_step;
   Box_II->Eig_Vec     = Ham_Info->Vector[0];
   Box_II->Eig_Val     = &Ham_Info->Value[0];
   Box_II->p_threads   = p_threads;
   strcpy(Box_II->Type, "CG");
   
   INVERSE_ITERATION(Box_II);
   
   free(Box_II);
   
   if (Ham_Info->diag_num >= 2) {
      printf("Underconstruction\n");
      exit(1);
   }
   
   Ham_Info->Error[0] = CONVERGE_CHECK(Ham_Info->Ham, Ham_Info->Vector[0], Ham_Info->Value[0], p_threads);
   Ham_Info->mem_ham  = (double)Ham_Info->Ham->Row[Ham_Info->Ham->row_dim]*12*pow(10,-9) + (double)Ham_Info->Ham->row_dim*8*pow(10,-9);
   Ham_Info->dim      = Ham_Info->Ham->row_dim;
   Time->inv_iter     = omp_get_wtime() - Time->inv_iter;
   
}
