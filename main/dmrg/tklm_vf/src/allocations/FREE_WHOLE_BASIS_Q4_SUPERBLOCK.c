//
//  FREE_WHOLE_BASIS_Q4_SUPERBLOCK.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/22.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q4_SUPERBLOCK(DMRG_WHOLE_BASIS_Q4 *Basis, MODEL_1DTKLM_VF *Model) {
   
   int sz_sign,tot_sz,del_Ne_1,del_Ne_2;
   int spin   = Model->spin_loc;
   int max_sz = spin*Model->tot_site + Model->tot_ele_1 + Model->tot_ele_2 + 4*spin + 2;
   
   for (del_Ne_1 = 0; del_Ne_1 <= 2; del_Ne_1++) {
      for (del_Ne_2 = 0; del_Ne_2 <= 2; del_Ne_2++) {
         for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
            for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
               free(Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               free(Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               free(Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               free(Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
            }
            free(Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]);
            free(Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]);
            free(Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]);
            free(Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]);
         }
         free(Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2]);
         free(Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2]);
         free(Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2]);
         free(Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2]);
      }
      free(Basis->LL_LLLRRRRL[del_Ne_1]);
      free(Basis->LR_LLLRRRRL[del_Ne_1]);
      free(Basis->RL_LLLRRRRL[del_Ne_1]);
      free(Basis->RR_LLLRRRRL[del_Ne_1]);
   }
   free(Basis->LL_LLLRRRRL);
   free(Basis->LR_LLLRRRRL);
   free(Basis->RL_LLLRRRRL);
   free(Basis->RR_LLLRRRRL);
   
   FREE_ARRAY_INT4(Basis->Dim, 3, 3, 2);
   
   free(Basis);
}
