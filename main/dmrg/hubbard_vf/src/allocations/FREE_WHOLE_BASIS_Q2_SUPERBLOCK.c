//
//  FREE_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DHUBBARD_VF *Model) {
   
   int i,j;
   
   int max_sz = Model->tot_site + Model->tot_ele + 2;
   
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
