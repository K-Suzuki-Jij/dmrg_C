#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LL_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M, double coeef, int *elem_num) {
   
   long i;
   
   int LL         = Basis->LL;
   int LR         = Basis->LR;
   int RR         = Basis->RR;
   int RL         = Basis->RL;
   int dim_RR     = Basis->dim_RR;
   int dim_onsite = Basis->dim_onsite;
   int t_elem_num = *elem_num;
   int inv_check,col_LL,inv;
   long inv_sup;
   double val;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   for (i = M->Row[LL]; i < M->Row[LL + 1]; i++) {
      col_LL    = M->Col[i];
      inv_sup   = (long)col_LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      inv       = Inv_LLLRRRRL[inv_sup];
      inv_check = A_Basis->Inv_LLLRRRRL[inv];
      val = coeef*M->Val[i];
      if (fabs(val) > zero && inv >= 0) {
         if (inv_check == -1) {
            A_Basis->LL_LLLRRRRL[t_elem_num]  = col_LL;
            A_Basis->LR_LLLRRRRL[t_elem_num]  = LR;
            A_Basis->RR_LLLRRRRL[t_elem_num]  = RR;
            A_Basis->RL_LLLRRRRL[t_elem_num]  = RL;
            A_Basis->Val_LLLRRRRL[t_elem_num] = val;
            A_Basis->Inv_LLLRRRRL[inv]        = t_elem_num;
            t_elem_num++;
         }
         else {
            A_Basis->Val_LLLRRRRL[inv_check] = A_Basis->Val_LLLRRRRL[inv_check] + val;
         }
      }
   }
   
   *elem_num = t_elem_num;
   
}
