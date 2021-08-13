//
//  FREE_BLOCK_MATRIX.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   
   FREE_CRS2(Block->Ham      , tot_site);
   FREE_CRS2(Block->SpL_RE   , tot_site);
   FREE_CRS2(Block->SmL_RE   , tot_site);
   FREE_CRS2(Block->SzL_RE   , tot_site);
   FREE_CRS2(Block->Even_RE  , tot_site);
   FREE_CRS2(Block->Odd_RE   , tot_site);
   FREE_CRS2(Block->Even_D_RE, tot_site);
   FREE_CRS2(Block->Odd_D_RE , tot_site);
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      FREE_CRS2(Block->SpL_LE   , tot_site);
      FREE_CRS2(Block->SmL_LE   , tot_site);
      FREE_CRS2(Block->SzL_LE   , tot_site);
      FREE_CRS2(Block->Even_LE  , tot_site);
      FREE_CRS2(Block->Odd_LE   , tot_site);
      FREE_CRS2(Block->Even_D_LE, tot_site);
      FREE_CRS2(Block->Odd_D_LE , tot_site);
   }
   
   FREE_CRS1(Block->SzL_On   );
   FREE_CRS1(Block->SpL_On   );
   FREE_CRS1(Block->SmL_On   );
   FREE_CRS1(Block->Even_On  );
   FREE_CRS1(Block->Odd_On   );
   FREE_CRS1(Block->Even_D_On);
   FREE_CRS1(Block->Odd_D_On );
   FREE_CRS1(Block->Ham_On   );
   
   FREE_ARRAY_INT2(Block->Tot_Parity, tot_site);
   FREE_ARRAY_INT2(Block->Tot_Ele, tot_site);
   FREE_ARRAY_INT1(Block->Dim);
   
   FREE_CRS1(Block->SxL_On );
   FREE_CRS1(Block->SzC_On );
   FREE_CRS1(Block->SxC_On );
   FREE_CRS1(Block->SCSL_On);
   FREE_CRS1(Block->NC_On  );
   
   free(Block->TM);
   free(Block->TM_D);
   
   free(Block->Basis_LL_LLLR);
   free(Block->Basis_LR_LLLR);
   
   FREE_ARRAY_INT3(Block->Basis_Inv_LLLR, tot_site/2, max_dim);
   FREE_ARRAY_INT1(Block->Dim_LLLR);
   
   free(Block);
   
}
