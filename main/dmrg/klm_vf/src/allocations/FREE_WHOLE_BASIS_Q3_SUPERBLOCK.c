//
//  FREE_WHOLE_BASIS_Q3_SUPERBLOCK.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/08/14.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q3_SUPERBLOCK(DMRG_WHOLE_BASIS_Q3 *Basis, MODEL_1DKLM_VF *Model) {
   
   int sz_sign,tot_sz,del_Ne;
   int max_sz     = Model->spin_loc*Model->tot_site + Model->tot_ele + 4*Model->spin_loc + 2;
   
   for (del_Ne = 0; del_Ne <= 2; del_Ne++) {
      for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
         for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
            free(Basis->LL_LLLRRRRL[del_Ne][sz_sign][tot_sz]);
            free(Basis->LR_LLLRRRRL[del_Ne][sz_sign][tot_sz]);
            free(Basis->RL_LLLRRRRL[del_Ne][sz_sign][tot_sz]);
            free(Basis->RR_LLLRRRRL[del_Ne][sz_sign][tot_sz]);
         }
         free(Basis->LL_LLLRRRRL[del_Ne][sz_sign]);
         free(Basis->LR_LLLRRRRL[del_Ne][sz_sign]);
         free(Basis->RL_LLLRRRRL[del_Ne][sz_sign]);
         free(Basis->RR_LLLRRRRL[del_Ne][sz_sign]);
      }
      free(Basis->LL_LLLRRRRL[del_Ne]);
      free(Basis->LR_LLLRRRRL[del_Ne]);
      free(Basis->RL_LLLRRRRL[del_Ne]);
      free(Basis->RR_LLLRRRRL[del_Ne]);
   }
   
   free(Basis->LL_LLLRRRRL);
   free(Basis->LR_LLLRRRRL);
   free(Basis->RL_LLLRRRRL);
   free(Basis->RR_LLLRRRRL);
   
   FREE_ARRAY_INT3(Basis->Dim, 3, 2);
   
   free(Basis);
}
