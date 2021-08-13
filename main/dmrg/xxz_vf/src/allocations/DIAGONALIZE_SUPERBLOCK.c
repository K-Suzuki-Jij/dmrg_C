//
//  DIAGONALIZE_SUPERBLOCK.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/27.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void DIAGONALIZE_SUPERBLOCK(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_ham = omp_get_wtime();
   
   CRS1 *Ham_LLLR      = GET_HAM_LLLR(System, Enviro, Dmrg_Basis, Dmrg_Status, Model);
   CRS1 *Ham_LRRL      = GET_HAM_LRRL(System, Enviro, Dmrg_Basis, Dmrg_Status, Model);
   CRS1 *Ham_RRRL      = GET_HAM_RRRL(System, Enviro, Dmrg_Basis, Dmrg_Status, Model);
   CRS1 *Ham_LRRL_Sign = GET_CRS1(Dmrg_Basis->dim_LLLR, 1);
   int *Ele_RR = GET_ARRAY_INT1(Dmrg_Basis->dim_RRRL);
   double eig_val = 0;
         
   Dmrg_Time->make_ham = omp_get_wtime() - Dmrg_Time->make_ham;
   
   Dmrg_Time->diag = omp_get_wtime();
   
   DMRG_BOX_LAN *Box_Lan = malloc(sizeof(DMRG_BOX_LAN));
   Box_Lan->Eig_Vec      = Dmrg_Status->GS_Vec;
   Box_Lan->eig_val      = &eig_val;
   Box_Lan->acc          = Dmrg_Param->diag_acc;
   Box_Lan->min_step     = Dmrg_Param->diag_min_step;
   Box_Lan->max_step     = Dmrg_Param->diag_max_step;
   Box_Lan->p_threads    = Model->p_threads;
   Box_Lan->M_LLLR       = Ham_LLLR;
   Box_Lan->M_LRRL       = Ham_LRRL;
   Box_Lan->M_LRRL_Sign  = Ham_LRRL_Sign;
   Box_Lan->M_RRRL       = Ham_RRRL;
   Box_Lan->Ele_RR       = Ele_RR;
   Box_Lan->Dmrg_Basis   = Dmrg_Basis;
   strcpy(Box_Lan->Guess, Dmrg_Param->Initial_Guess);
   strcpy(Box_Lan->Type , "Normal");

   DMRG_LANCZOS_SLOW_OBC(Box_Lan);
   
   Dmrg_Time->diag = omp_get_wtime() - Dmrg_Time->diag;
   
   Dmrg_Time->inv_iter = omp_get_wtime();

   DMRG_BOX_II *Box_II = malloc(sizeof(DMRG_BOX_II));
   Box_II->ii_acc      = Dmrg_Param->inv_iter_acc;
   Box_II->ii_diag_add = Dmrg_Param->inv_iter_diag_add;
   Box_II->ii_max_step = Dmrg_Param->inv_iter_max_step;
   Box_II->cg_acc      = Dmrg_Param->cg_acc;
   Box_II->cg_max_step = Dmrg_Param->cg_max_step;
   Box_II->Eig_Vec     = Box_Lan->Eig_Vec;
   Box_II->eig_val     = Box_Lan->eig_val;
   Box_II->p_threads   = Model->p_threads;
   Box_II->M_LLLR      = Ham_LLLR;
   Box_II->M_LRRL      = Ham_LRRL;
   Box_II->M_LRRL_Sign = Ham_LRRL_Sign;
   Box_II->M_RRRL      = Ham_RRRL;
   Box_II->Ele_RR      = Ele_RR;
   Box_II->Dmrg_Basis  = Dmrg_Basis;
   
   DMRG_INVERSE_ITERATION_OBC(Box_II);

   Dmrg_Status->gs_val   = eig_val;
   Dmrg_Status->gs_error = Box_II->error;

   Dmrg_Time->inv_iter = omp_get_wtime() - Dmrg_Time->inv_iter;

   FREE_CRS1(Ham_LLLR);
   FREE_CRS1(Ham_LRRL);
   FREE_CRS1(Ham_RRRL);
   FREE_CRS1(Ham_LRRL_Sign);
   FREE_ARRAY_INT1(Ele_RR);
   free(Box_Lan);
   free(Box_II);
   
}
