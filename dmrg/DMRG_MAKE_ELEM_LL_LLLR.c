#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LL_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR,  CRS1 *M, double coeef, int *elem_num) {
   
   long i;
   
   int LL = Basis->LL;
   int LR = Basis->LR;
   int t_elem_num = *elem_num;
   
   int inv_check,col_LL,inv;
   double val;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }

   for (i = M->Row[LL]; i < M->Row[LL + 1]; i++) {
      col_LL    = M->Col[i];
      inv       = Inv_LLLR[col_LL][LR];
      inv_check = A_Basis->Inv_LLLRRRRL[inv];
      val       = coeef*M->Val[i];
      if (fabs(val) > zero && inv >= 0) {
         if (inv_check == -1) {
            A_Basis->LL_LLLRRRRL[t_elem_num]   = col_LL;
            A_Basis->LR_LLLRRRRL[t_elem_num]   = LR;
            A_Basis->Val_LLLRRRRL[t_elem_num]  = val;
            A_Basis->Inv_LLLRRRRL[inv]         = t_elem_num;
            t_elem_num++;
         }
         else {
            A_Basis->Val_LLLRRRRL[inv_check] = A_Basis->Val_LLLRRRRL[inv_check] + val;
         }
      }
   }
   
   *elem_num = t_elem_num;
   
}
