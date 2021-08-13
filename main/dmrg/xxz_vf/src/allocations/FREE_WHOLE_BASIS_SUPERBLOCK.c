//
//  FREE_WHOLE_BASIS_SUPERBLOCK.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/06.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DXXZ_VF *Model) {
   
   int sz_sign,tot_sz;
   int max_sz = Model->spin*Model->tot_site;
   
   for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
      for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
         free(Basis->LL_LLLRRRRL[sz_sign][tot_sz]);
         free(Basis->LR_LLLRRRRL[sz_sign][tot_sz]);
         free(Basis->RL_LLLRRRRL[sz_sign][tot_sz]);
         free(Basis->RR_LLLRRRRL[sz_sign][tot_sz]);
      }
      free(Basis->LL_LLLRRRRL[sz_sign]);
      free(Basis->LR_LLLRRRRL[sz_sign]);
      free(Basis->RL_LLLRRRRL[sz_sign]);
      free(Basis->RR_LLLRRRRL[sz_sign]);
   }
   
   free(Basis->LL_LLLRRRRL);
   free(Basis->LR_LLLRRRRL);
   free(Basis->RL_LLLRRRRL);
   free(Basis->RR_LLLRRRRL);
   
   FREE_ARRAY_INT2(Basis->Dim, 2);
   
   free(Basis);
}
