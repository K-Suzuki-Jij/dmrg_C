#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LLRR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LL, CRS1 *M_RR, double coeef, int *Ele_LL, int *Ele_LR, int *elem_num, char Type[], char Sign_Flag[]) {
   
   long i,j;
   
   int LL         = Basis->LL;
   int LR         = Basis->LR;
   int RR         = Basis->RR;
   int RL         = Basis->RL;
   int dim_RR     = Basis->dim_RR;
   int dim_onsite = Basis->dim_onsite;
   int t_elem_num = *elem_num;
   int inv_check,sign,col_LL,col_RR,inv;
   long inv_sup;
   double val,val_LL,val_RR;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   if (strcmp(Sign_Flag, "Yes") == 0) {
      if (strcmp(Type, "LL_RR") == 0) {
         if ((Ele_LL[LL] + Ele_LR[LR])%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (strcmp(Type, "RR_LL") == 0) {
         if ((Ele_LL[LL] + Ele_LR[LR])%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         printf("Error in DMRG_MAKE_ELEM_LLRR_LLLRRRRL\n");
         exit(1);
      }
   }
   else {
      sign = 1;
   }
   
   for (i = M_LL->Row[LL]; i < M_LL->Row[LL + 1]; i++) {
      col_LL = M_LL->Col[i];
      val_LL = M_LL->Val[i];
      for (j = M_RR->Row[RR]; j < M_RR->Row[RR + 1]; j++) {
         col_RR    = M_RR->Col[j];
         val_RR    = M_RR->Val[j];
         inv_sup   = (long)col_LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + col_RR*dim_onsite + RL;
         inv       = Inv_LLLRRRRL[inv_sup];
         inv_check = A_Basis->Inv_LLLRRRRL[inv];
         val       = sign*coeef*val_LL*val_RR;
         if (fabs(val) > zero && inv >= 0) {
            if (inv_check == -1) {
               A_Basis->LL_LLLRRRRL[t_elem_num]  = col_LL;
               A_Basis->LR_LLLRRRRL[t_elem_num]  = LR;
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
   }
   
   *elem_num = t_elem_num;
   
}
