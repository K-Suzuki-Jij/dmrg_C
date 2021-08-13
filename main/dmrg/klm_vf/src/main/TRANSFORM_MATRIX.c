//
//  TRANSFORM_MATRIX.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/07.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void TRANSFORM_MATRIX(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DKLM_VF *Model, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->trans_main = omp_get_wtime();
   
   int LL_site  = Dmrg_Status->LL_site;
   int dim_LLLR = Dmrg_Basis->dim_LLLR;
   int elem_num = dim_LLLR*Dmrg_System->dim_renorm;
   CCS1 *T_M    = Dmrg_System->Trans_Matrix;
   CRS1 *T_MD   = Dmrg_System->Trans_Matrix_Dagger;
   CRS1 *M_CRS  = GET_CRS1(dim_LLLR, elem_num);
   CCS1 *M_CCS  = GET_CCS1(dim_LLLR, elem_num);
   double *V    = GET_ARRAY_DOUBLE1(dim_LLLR);
   
   CRS1 *Ham_LLLR = GET_HAM_LLLR(Block_System, Block_Enviro, Dmrg_Basis, Dmrg_Status, Model);
   CRS_CRS_CCS_PRODUCT(T_MD, Ham_LLLR, T_M, Block_System->Ham[LL_site + 1], M_CCS, V);
   FREE_CRS1(Ham_LLLR);
   
   DMRG_BASIS_LLLR *Dmrg_Basis_LLLR = malloc(sizeof(DMRG_BASIS_LLLR));
   Dmrg_Basis_LLLR->dim_LLLR = Dmrg_Basis->dim_LLLR;
   Dmrg_Basis_LLLR->LL_LLLR  = Dmrg_Basis->LL_LLLR;
   Dmrg_Basis_LLLR->LR_LLLR  = Dmrg_Basis->LR_LLLR;
   Dmrg_Basis_LLLR->Inv_LLLR = Dmrg_Basis->Inv_LLLR;
   
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->CUp_On    , Block_System->CUp_RE[LL_site + 1]    , Block_System->Tot_Ele[LL_site], "LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->CUp_D_On  , Block_System->CUp_D_RE[LL_site + 1]  , Block_System->Tot_Ele[LL_site], "LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->CDown_On  , Block_System->CDown_RE[LL_site + 1]  , Block_System->Tot_Ele[LL_site], "LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->CDown_D_On, Block_System->CDown_D_RE[LL_site + 1], Block_System->Tot_Ele[LL_site], "LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->SpL_On    , Block_System->SpL_RE[LL_site + 1]    , NULL                          , "LR", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->SmL_On    , Block_System->SmL_RE[LL_site + 1]    , NULL                          , "LR", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   DMRG_EXTEND_AND_REDUCE(NULL, Block_System->SzL_On    , Block_System->SzL_RE[LL_site + 1]    , NULL                          , "LR", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      DMRG_EXTEND_AND_REDUCE(Block_System->CUp_LE[LL_site]    , NULL, Block_System->CUp_LE[LL_site + 1]    , Block_System->Tot_Ele[LL_site], "LL", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->CUp_D_LE[LL_site]  , NULL, Block_System->CUp_D_LE[LL_site + 1]  , Block_System->Tot_Ele[LL_site], "LL", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->CDown_LE[LL_site]  , NULL, Block_System->CDown_LE[LL_site + 1]  , Block_System->Tot_Ele[LL_site], "LL", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->CDown_D_LE[LL_site], NULL, Block_System->CDown_D_LE[LL_site + 1], Block_System->Tot_Ele[LL_site], "LL", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->SpL_LE[LL_site]    , NULL, Block_System->SpL_LE[LL_site + 1]    , NULL                          , "LL", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->SmL_LE[LL_site]    , NULL, Block_System->SmL_LE[LL_site + 1]    , NULL                          , "LL", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
      DMRG_EXTEND_AND_REDUCE(Block_System->SzL_LE[LL_site]    , NULL, Block_System->SzL_LE[LL_site + 1]    , NULL                          , "LL", "No" , T_M, T_MD, M_CRS, M_CCS, V, Dmrg_Basis_LLLR);
   }
   
   free(Dmrg_Basis_LLLR);
   FREE_CRS1(M_CRS);
   FREE_CCS1(M_CCS);
   FREE_ARRAY_DOUBLE1(V);
   
   int basis;
   
   Block_System->Dim[LL_site + 1] = Dmrg_System->dim_renorm;
   
   for (basis = 0; basis < Dmrg_System->dim_renorm; basis++) {
      Block_System->Tot_Ele[LL_site + 1][basis] = Dmrg_System->Q_Number1[basis];
      Block_System->Tot_Sz[LL_site + 1][basis]  = Dmrg_System->Q_Number2[basis];
   }
   
   Dmrg_Status->percent_LL = FIND_PERCENTAGE_LL(Block_System, Model);
   
   Dmrg_Time->trans_main = omp_get_wtime() - Dmrg_Time->trans_main;
   
}
