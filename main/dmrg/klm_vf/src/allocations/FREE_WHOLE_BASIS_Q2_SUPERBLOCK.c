//
//  FREE_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/07.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DKLM_VF *Model) {
   
   int i,j;
   
   int max_sz = Model->spin_loc*Model->tot_site + Model->tot_ele + 4*Model->spin_loc + 2;

   
   for (i = 0; i < 2; i++) {
      for (j = 0; j < max_sz + 1; j++) {
         free(Basis->LL_LLLRRRRL[i][j]);
         free(Basis->LR_LLLRRRRL[i][j]);
         free(Basis->RL_LLLRRRRL[i][j]);
         free(Basis->RR_LLLRRRRL[i][j]);
      }
      free(Basis->LL_LLLRRRRL[i]);
      free(Basis->LR_LLLRRRRL[i]);
      free(Basis->RL_LLLRRRRL[i]);
      free(Basis->RR_LLLRRRRL[i]);
   }
   free(Basis->LL_LLLRRRRL);
   free(Basis->LR_LLLRRRRL);
   free(Basis->RL_LLLRRRRL);
   free(Basis->RR_LLLRRRRL);
   
   FREE_ARRAY_INT2(Basis->Dim, 2);
   
   free(Basis);
}
