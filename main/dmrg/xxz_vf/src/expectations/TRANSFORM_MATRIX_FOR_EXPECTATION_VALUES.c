//
//  TRANSFORM_MATRIX_FOR_EXPECTATION_VALUES.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/21.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void TRANSFORM_MATRIX_FOR_EXPECTATION_VALUES(BLOCK *Block_System, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   int c1 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 3);
   int c3 = (Dmrg_Status->LL_site >= Model->tot_site/2 - 2);
   int c4 = (strcmp(Dmrg_Status->Enviro_Copy, "yes") == 0);
   int c5 = (strcmp(Dmrg_Status->Enviro_Copy, "yes") != 0);
   
   if (c3 || (c4 && c1) || (c5 && c2)) {
      Dmrg_Time->trans_exp = 0;
      Dmrg_Status->percent_LL = FIND_PERCENTAGE_LL(Block_System, Model);
      return;
   }
   
   Dmrg_Time->trans_exp = omp_get_wtime();
   
   int LL_site  = Dmrg_Status->LL_site;
   int dim_LLLR = Dmrg_Basis->dim_LLLR;
   int elem_num = dim_LLLR*Dmrg_System->dim_renorm;
   CCS1 *T_M    = Dmrg_System->Trans_Matrix;
   CRS1 *T_MD   = Dmrg_System->Trans_Matrix_Dagger;
   CRS1 **M_CRS = GET_CRS2(Model->p_threads, dim_LLLR, elem_num);
   CCS1 **M_CCS = GET_CCS2(Model->p_threads, dim_LLLR, elem_num);
   double **V   = GET_ARRAY_DOUBLE2(Model->p_threads, dim_LLLR);
   
   int r,site = 0,thread_num;
   
   

   if (LL_site == 0) {
      ONSITE_SZ_SZBASIS_HB(Model->spin, Block_System->Sz[0], 1.0);
      ONSITE_SX_SZBASIS_HB(Model->spin, Block_System->Sx[0], 1.0);
      ONSITE_SZSZ_SZBASIS_HB(Model->spin, Block_System->SzSz[0], 1.0);
      ONSITE_SXSX_SZBASIS_HB(Model->spin, Block_System->SxSx[0], 1.0);

      //Correlation Functions
      if (Model->cf_origin == 0) {
         ONSITE_SZSZ_SZBASIS_HB(Model->spin, Block_System->Sz_CF[0], 1.0);
         ONSITE_SXSX_SZBASIS_HB(Model->spin, Block_System->Sx_CF[0], 1.0);
      }
   }
   
   DMRG_BASIS_LLLR *Dmrg_Basis_LLLR = malloc(sizeof(DMRG_BASIS_LLLR));
   Dmrg_Basis_LLLR->dim_LLLR = Dmrg_Basis->dim_LLLR;
   Dmrg_Basis_LLLR->LL_LLLR = Dmrg_Basis->LL_LLLR;
   Dmrg_Basis_LLLR->LR_LLLR = Dmrg_Basis->LR_LLLR;
   Dmrg_Basis_LLLR->Inv_LLLR = Dmrg_Basis->Inv_LLLR;

   //Correlation Functions
   if (LL_site >= Model->cf_origin) {
      r = LL_site - Model->cf_origin;

      DMRG_EXTEND_AND_REDUCE(Block_System->Sz[Model->cf_origin], Block_System->Sz_On, Block_System->Sz_CF[r + 1], NULL, "LL_LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->Sx[Model->cf_origin], Block_System->Sx_On, Block_System->Sx_CF[r + 1], NULL, "LL_LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);

      if (r == 0) {
         DMRG_EXTEND_AND_REDUCE(Block_System->SzSz[LL_site], NULL, Block_System->Sz_CF[r], NULL, "LL", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
         DMRG_EXTEND_AND_REDUCE(Block_System->SxSx[LL_site], NULL, Block_System->Sx_CF[r], NULL, "LL", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
      }
      else {
#pragma omp parallel for private (thread_num) num_threads (Model->p_threads)
         for (r = 0; r <= LL_site - Model->cf_origin; r++){
            thread_num = omp_get_thread_num();
            DMRG_EXTEND_AND_REDUCE(Block_System->Sz_CF[r], NULL, Block_System->Sz_CF[r], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
            DMRG_EXTEND_AND_REDUCE(Block_System->Sx_CF[r], NULL, Block_System->Sx_CF[r], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
         }
      }
   }

#pragma omp parallel for private (thread_num) num_threads (Model->p_threads)
   for (site = 0; site <= LL_site; site++) {
      thread_num = omp_get_thread_num();
      DMRG_EXTEND_AND_REDUCE(Block_System->Sz[site]  , NULL, Block_System->Sz[site]  , NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->SzSz[site], NULL, Block_System->SzSz[site], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->Sx[site]  , NULL, Block_System->Sx[site]  , NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->SxSx[site], NULL, Block_System->SxSx[site], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Dmrg_Basis_LLLR);
   }
   
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->Sz_On  , Block_System->Sz[LL_site + 1]  , NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->SzSz_On, Block_System->SzSz[LL_site + 1], NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->Sx_On  , Block_System->Sx[LL_site + 1]  , NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->SxSx_On, Block_System->SxSx[LL_site + 1], NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Dmrg_Basis_LLLR);
   
   free(Dmrg_Basis_LLLR);
   FREE_CRS2(M_CRS, Model->p_threads);
   FREE_CCS2(M_CCS, Model->p_threads);
   FREE_ARRAY_DOUBLE2(V, Model->p_threads);
   
   Dmrg_Status->percent_LL = FIND_PERCENTAGE_LL(Block_System, Model);
   
   Dmrg_Time->trans_exp = omp_get_wtime() - Dmrg_Time->trans_exp;
   
}
