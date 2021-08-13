//
//  FREE_BLOCK_MATRIX.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   
   FREE_CRS2(Block->Ham         , tot_site);
   FREE_CRS2(Block->CUp_RE      , tot_site);
   FREE_CRS2(Block->CDown_RE    , tot_site);
   FREE_CRS2(Block->CUp_D_RE    , tot_site);
   FREE_CRS2(Block->CDown_D_RE  , tot_site);
   FREE_CRS2(Block->NC_Up_RE    , tot_site);
   FREE_CRS2(Block->NC_Down_RE  , tot_site);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      FREE_CRS2(Block->CUp_LE      , tot_site);
      FREE_CRS2(Block->CDown_LE    , tot_site);
      FREE_CRS2(Block->CUp_D_LE    , tot_site);
      FREE_CRS2(Block->CDown_D_LE  , tot_site);
      FREE_CRS2(Block->NC_Up_LE    , tot_site);
      FREE_CRS2(Block->NC_Down_LE  , tot_site);
   }
   
   FREE_CRS1(Block->CUp_On    );
   FREE_CRS1(Block->CDown_On  );
   FREE_CRS1(Block->CUp_D_On  );
   FREE_CRS1(Block->CDown_D_On);
   FREE_CRS1(Block->Ham_On    );
   FREE_CRS1(Block->NC_Up_On  );
   FREE_CRS1(Block->NC_Down_On);

   FREE_ARRAY_INT2(Block->Tot_Sz , tot_site);
   FREE_ARRAY_INT2(Block->Tot_Ele, tot_site);
   FREE_ARRAY_INT1(Block->Dim);
   
   FREE_CRS1(Block->SzC_On );
   FREE_CRS1(Block->SxC_On );
   FREE_CRS1(Block->NC_On  );
   FREE_CRS1(Block->DO_On  );


   
   free(Block->TM);
   free(Block->TM_D);
   
   free(Block->Basis_LL_LLLR);
   free(Block->Basis_LR_LLLR);
   
   FREE_ARRAY_INT3(Block->Basis_Inv_LLLR, tot_site/2, max_dim);
   FREE_ARRAY_INT1(Block->Dim_LLLR);
   
   
   free(Block);
   
}
