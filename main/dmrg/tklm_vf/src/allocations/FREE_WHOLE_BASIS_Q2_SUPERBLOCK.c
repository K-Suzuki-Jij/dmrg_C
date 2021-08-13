//
//  FREE_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/24.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis, MODEL_1DTKLM_VF *Model) {
   
   int sz_sign,tot_sz;
   int spin   = Model->spin_loc;
   int max_sz = spin*Model->tot_site + Model->tot_ele_1 + Model->tot_ele_2 + 4*spin + 2;
   
   
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
