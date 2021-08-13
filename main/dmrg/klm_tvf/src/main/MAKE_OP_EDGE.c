//
//  MAKE_OP_EDGE.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DKLM_TVF *Model) {
   
   int basis;
   for (basis = 0; basis < Model->dim_onsite; basis++) {
      Block->Tot_Parity[0][basis] = ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(basis, Model->spin_loc);
      Block->Tot_Ele[0][basis]    = ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(basis, Model->spin_loc);
   }
   
   Block->Dim[0] = Model->dim_onsite;
   
   COPY_CRS1(Block->Ham_On    , Block->Ham[0]       , 1);
   COPY_CRS1(Block->SpL_On    , Block->SpL_RE[0]    , 1);
   COPY_CRS1(Block->SmL_On    , Block->SmL_RE[0]    , 1);
   COPY_CRS1(Block->SzL_On    , Block->SzL_RE[0]    , 1);
   COPY_CRS1(Block->Even_On   , Block->Even_RE[0]   , 1);
   COPY_CRS1(Block->Odd_On    , Block->Odd_RE[0]    , 1);
   COPY_CRS1(Block->Even_D_On , Block->Even_D_RE[0] , 1);
   COPY_CRS1(Block->Odd_D_On  , Block->Odd_D_RE[0]  , 1);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      COPY_CRS1(Block->SpL_On    , Block->SpL_LE[0]    , 1);
      COPY_CRS1(Block->SmL_On    , Block->SmL_LE[0]    , 1);
      COPY_CRS1(Block->SzL_On    , Block->SzL_LE[0]    , 1);
      COPY_CRS1(Block->Even_On   , Block->Even_LE[0]   , 1);
      COPY_CRS1(Block->Odd_On    , Block->Odd_LE[0]    , 1);
      COPY_CRS1(Block->Even_D_On , Block->Even_D_LE[0] , 1);
      COPY_CRS1(Block->Odd_D_On  , Block->Odd_D_LE[0]  , 1);
   }
 
   
}
