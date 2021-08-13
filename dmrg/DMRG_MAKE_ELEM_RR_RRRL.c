#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_RR_RRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_RRRL, CRS1 *M, double coeef, int *elem_num) {
   
   long i;
   
   int RR  = Basis->RR;
   int RL  = Basis->RL;
   int t_elem_num = *elem_num;
   int inv_check,col_RR,inv;
   double val;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   for (i = M->Row[RR]; i < M->Row[RR + 1]; i++) {
      col_RR    = M->Col[i];
      inv       = Inv_RRRL[col_RR][RL];
      inv_check = A_Basis->Inv_LLLRRRRL[inv];
      val       = coeef*M->Val[i];
      if (fabs(val) > zero && inv >= 0) {
         if (inv_check == -1) {
            A_Basis->RR_LLLRRRRL[t_elem_num]  = col_RR;
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
