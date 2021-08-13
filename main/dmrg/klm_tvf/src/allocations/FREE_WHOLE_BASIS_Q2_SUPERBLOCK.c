//
//  FREE_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2_SUPERBLOCK(DMRG_WHOLE_BASIS_Q2 *Basis) {
   
   FREE_ARRAY_SINT3(Basis->LL_LLLRRRRL, 3, 2);
   FREE_ARRAY_SINT3(Basis->LR_LLLRRRRL, 3, 2);
   FREE_ARRAY_SINT3(Basis->RL_LLLRRRRL, 3, 2);
   FREE_ARRAY_SINT3(Basis->RR_LLLRRRRL, 3, 2);
   FREE_ARRAY_INT2(Basis->Dim, 3);
   
   free(Basis);
}
