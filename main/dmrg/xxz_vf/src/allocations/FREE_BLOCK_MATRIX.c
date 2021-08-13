//
//  FREE_BLOCK_MATRIX.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/24.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_BLOCK_MATRIX(BLOCK *Block, MODEL_1DXXZ_VF *Model) {
   
   int tot_site   = Model->tot_site;
   int cf_length  = Model->tot_site/2 - 1 - Model->cf_origin;
   
   FREE_CRS2(Block->Ham  , tot_site);
   FREE_CRS2(Block->Sz_RE, tot_site);
   FREE_CRS2(Block->Sp_RE, tot_site);
   FREE_CRS2(Block->Sm_RE, tot_site);
   
   if (strcmp(Model->BC, "PBC") == 0) {
      FREE_CRS2(Block->Sz_LE, tot_site);
      FREE_CRS2(Block->Sp_LE, tot_site);
      FREE_CRS2(Block->Sm_LE, tot_site);
   }
   
   FREE_CRS1(Block->Sz_On  );
   FREE_CRS1(Block->Sp_On  );
   FREE_CRS1(Block->Sm_On  );
   FREE_CRS1(Block->Ham_On );
   FREE_CRS1(Block->Sx_On  );
   FREE_CRS1(Block->SzSz_On);
   FREE_CRS1(Block->SxSx_On);
   
   FREE_CRS2(Block->Sz  , tot_site/2);
   FREE_CRS2(Block->SzSz, tot_site/2);
   FREE_CRS2(Block->Sx  , tot_site/2);
   FREE_CRS2(Block->SxSx, tot_site/2);
   
   FREE_CRS2(Block->Sz_CF, cf_length);
   FREE_CRS2(Block->Sx_CF, cf_length);
   
   FREE_ARRAY_INT2(Block->Tot_Sz, tot_site);
   FREE_ARRAY_INT1(Block->Dim);
   
   free(Block);
   
   
}
