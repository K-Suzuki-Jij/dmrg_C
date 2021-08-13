//
//  FREE_WHOLE_BASIS_Q1_SUPERBLOCK.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/17.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q1_SUPERBLOCK(DMRG_WHOLE_BASIS_Q1 *Basis) {
   
   FREE_ARRAY_SINT2(Basis->LL_LLLRRRRL, 2);
   FREE_ARRAY_SINT2(Basis->LR_LLLRRRRL, 2);
   FREE_ARRAY_SINT2(Basis->RL_LLLRRRRL, 2);
   FREE_ARRAY_SINT2(Basis->RR_LLLRRRRL, 2);
   FREE_ARRAY_INT1(Basis->Dim);
   
   free(Basis);
}
