#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LR_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR, CRS1 *M, double coeef, int *elem_num) {
   
   long i;
   
   int LL = Basis->LL;
   int LR = Basis->LR;
   int t_elem_num = *elem_num;
   
   int inv_check,col_LR,inv;
   double val;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   for (i = M->Row[LR]; i < M->Row[LR + 1]; i++) {
      col_LR    = M->Col[i];
      inv       = Inv_LLLR[LL][col_LR];
      inv_check = A_Basis->Inv_LLLRRRRL[inv];
      val       = coeef*M->Val[i];
      if (fabs(val) > zero && inv >= 0) {
         if (inv_check == -1) {
            A_Basis->LL_LLLRRRRL[t_elem_num]  = LL;
            A_Basis->LR_LLLRRRRL[t_elem_num]  = col_LR;
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
