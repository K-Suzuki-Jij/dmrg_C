#include "dmrg.h"

void DMRG_MAKE_ELEM_ZERO_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, int *elem_num) {
   
   int LL         = Basis->LL;
   int LR         = Basis->LR;
   int RR         = Basis->RR;
   int RL         = Basis->RL;
   int dim_RR     = Basis->dim_RR;
   int dim_onsite = Basis->dim_onsite;
   int t_elem_num = *elem_num;
   int inv_check,inv;
   long inv_sup;

   inv_sup   = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
   inv = Inv_LLLRRRRL[inv_sup];
   inv_check = A_Basis->Inv_LLLRRRRL[inv];
   if (inv_check == -1 && inv >= 0) {
      A_Basis->LL_LLLRRRRL[t_elem_num]  = LL;
      A_Basis->LR_LLLRRRRL[t_elem_num]  = LR;
      A_Basis->RR_LLLRRRRL[t_elem_num]  = RR;
      A_Basis->RL_LLLRRRRL[t_elem_num]  = RL;
      A_Basis->Val_LLLRRRRL[t_elem_num] = 0.0;
      A_Basis->Inv_LLLRRRRL[inv]        = t_elem_num;
      t_elem_num++;
   }
   
   *elem_num = t_elem_num;
   
}

