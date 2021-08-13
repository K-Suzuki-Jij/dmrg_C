//
//  MAKE_OP_EDGE.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DHUBBARD_VF *Model) {
   
   int basis;
   for (basis = 0; basis < Model->dim_onsite; basis++) {
      Block->Tot_Sz[0][basis]  = ONSITE_FIND_SITE_SZ_SZBASIS_HUBBARD(basis);
      Block->Tot_Ele[0][basis] = ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD(basis);
   }
   
   Block->Dim[0] = Model->dim_onsite;
   
   COPY_CRS1(Block->Ham_On    , Block->Ham[0]         , 1);
   COPY_CRS1(Block->CUp_On    , Block->CUp_RE[0]      , 1);
   COPY_CRS1(Block->CDown_On  , Block->CDown_RE[0]    , 1);
   COPY_CRS1(Block->CUp_D_On  , Block->CUp_D_RE[0]    , 1);
   COPY_CRS1(Block->CDown_D_On, Block->CDown_D_RE[0]  , 1);
   COPY_CRS1(Block->NC_Up_On  , Block->NC_Up_RE[0]    , 1);
   COPY_CRS1(Block->NC_Down_On, Block->NC_Down_RE[0]  , 1);

   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      COPY_CRS1(Block->CUp_On    , Block->CUp_LE[0]      , 1);
      COPY_CRS1(Block->CDown_On  , Block->CDown_LE[0]    , 1);
      COPY_CRS1(Block->CUp_D_On  , Block->CUp_D_LE[0]    , 1);
      COPY_CRS1(Block->CDown_D_On, Block->CDown_D_LE[0]  , 1);
      COPY_CRS1(Block->NC_Up_On  , Block->NC_Up_LE[0]    , 1);
      COPY_CRS1(Block->NC_Down_On, Block->NC_Down_LE[0]  , 1);
   }
   
   //For SSD
   MATRIX_CONSTAN_MULTIPLICATION_CRS1(Block->Ham[0], DMRG_SSD_COEFF(0, Model->tot_site, Model->BC, "LL", "Onsite"), 1);
   
}
