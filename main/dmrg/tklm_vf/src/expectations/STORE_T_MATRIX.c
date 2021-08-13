//
//  STORE_T_MATRIX.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void STORE_T_MATRIX(DMRG_BASIS *Basis, DMRG_SYSTEM_INFO *Dmrg_System, BLOCK *System, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, int LL_site, int dim_onsite) {
   
   int c1 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 3);
   int c3 = (Dmrg_Status->LL_site >= Model->tot_site/2 - 2);
   int c4 = (strcmp(Dmrg_Status->Enviro_Copy, "Yes") == 0);
   int c5 = (strcmp(Dmrg_Status->Enviro_Copy, "No") == 0);
   
   if (c3 || (c4 && c1) || (c5 && c2)) {
      return;
   }
   
   int dim       = Dmrg_System->Trans_Matrix_Dagger->row_dim;
   long elem_num = Dmrg_System->Trans_Matrix_Dagger->Row[dim];
   int dim_LLLR  = Basis->dim_LLLR;
   
   System->TM_D[LL_site] = GET_CRS1(dim, elem_num);
   System->TM[LL_site]   = GET_CCS1(dim, elem_num);
   
   COPY_CRS1(Dmrg_System->Trans_Matrix_Dagger, System->TM_D[LL_site], 1);
   COPY_CCS1(Dmrg_System->Trans_Matrix       , System->TM[LL_site]  , 1);
   
   System->Dim_LLLR[LL_site] = dim_LLLR;
   System->Basis_LL_LLLR[LL_site] = GET_ARRAY_SINT1(dim_LLLR);
   System->Basis_LR_LLLR[LL_site] = GET_ARRAY_SINT1(dim_LLLR);
   
   int i,j;
   for (i = 0; i < dim_LLLR; i++) {
      System->Basis_LL_LLLR[LL_site][i] = Basis->LL_LLLR[i];
      System->Basis_LR_LLLR[LL_site][i] = Basis->LR_LLLR[i];
   }
   
   for (i = 0; i < System->Dim[LL_site]; i++) {
      for (j = 0; j < dim_onsite; j++) {
         System->Basis_Inv_LLLR[LL_site][i][j] = Basis->Inv_LLLR[i][j];
      }
   }
   
}
