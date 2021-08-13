//
//  MAKE_OP_EDGE.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/18.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_OP_EDGE(BLOCK *Block, MODEL_1DXXZ_VF *Model) {
   
   int basis;
   for (basis = 0; basis < Model->dim_onsite; basis++) {
      Block->Tot_Sz[0][basis] = FIND_SITE_SZ(basis, Model->spin);
   }
   
   Block->Dim[0] = Model->dim_onsite;
   
   COPY_CRS1(Block->Ham_On, Block->Ham[0]   , 1);
   
   COPY_CRS1(Block->Sz_On , Block->Sz_RE[0] , 1);
   COPY_CRS1(Block->Sp_On , Block->Sp_RE[0] , 1);
   COPY_CRS1(Block->Sm_On , Block->Sm_RE[0] , 1);
   
   if (strcmp(Model->BC, "PBC") == 0) {
      COPY_CRS1(Block->Sz_On, Block->Sz_LE[0], 1);
      COPY_CRS1(Block->Sp_On, Block->Sp_LE[0], 1);
      COPY_CRS1(Block->Sm_On, Block->Sm_LE[0], 1);
   }
   
}
