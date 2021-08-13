#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "dmrg.h"

void DMRG_DIAGONALIZE_SUPERBLOCK_SYM(CRS1 *Ham, DMRG_TIME *Time, DMRG_PARAMETER *Param, DMRG_STATUS *Dmrg_Status, int p_threads) {
   
   //Diagonalize Hamiltonian
   Time->diag = omp_get_wtime();
   
   double eig_val = 0;
   
   BOX_LAN *Box_Lan        = malloc(sizeof(BOX_LAN));
   Box_Lan->M              = Ham;
   Box_Lan->Eig_Vec        = Dmrg_Status->GS_Vec;
   Box_Lan->eig_val        = &eig_val;
   Box_Lan->acc            = Param->diag_acc;
   Box_Lan->min_step       = Param->diag_min_step;
   Box_Lan->max_step       = Param->diag_max_step;
   Box_Lan->max_block_step = Param->diag_max_step;
   Box_Lan->p_threads      = p_threads;
   strcpy(Box_Lan->Type, Param->Lan_Con);
   strcpy(Box_Lan->Guess,"No");
   
   if (strcmp(Param->Diag_Method, "Lanczos") == 0) {
      LANCZOS_SYM(Box_Lan);
   }
   else if (strcmp(Param->Diag_Method, "Lanczos_Slow") == 0) {
      LANCZOS_SLOW_SYM(Box_Lan);
   }
   else {
      printf("Error in DMRG_DIAGONALIZE_HAMILTONIAN\n");
      exit(1);
   }
   
   Time->diag = omp_get_wtime() - Time->diag;
   
   //Inverse Iteration
   Time->inv_iter = omp_get_wtime();
   
   BOX_II *Box_II      = malloc(sizeof(BOX_II));
   Box_II->M           = Ham;
   Box_II->ii_acc      = Param->inv_iter_acc;
   Box_II->ii_diag_add = Param->inv_iter_diag_add;
   Box_II->ii_max_step = Param->inv_iter_max_step;
   Box_II->cg_acc      = Param->cg_acc;
   Box_II->cg_max_step = Param->cg_max_step;
   Box_II->Eig_Vec     = Dmrg_Status->GS_Vec;
   Box_II->Eig_Val     = &eig_val;
   Box_II->p_threads   = p_threads;
   strcpy(Box_II->Type, Param->II_Type);
   
   INVERSE_ITERATION_SYM(Box_II);

   Dmrg_Status->gs_val   = eig_val;
   Dmrg_Status->gs_error = Box_II->error;
   free(Box_II);
   free(Box_Lan);

   Time->inv_iter = omp_get_wtime() - Time->inv_iter;

}
