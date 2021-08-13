#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LRRR_LLLRRRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int *Inv_LLLRRRRL, CRS1 *M_LR, CRS1 *M_RR, double coeef, int *Ele_LR, int *elem_num, char Type[], char Sign_Flag[]) {
   
   long i,j;
   
   int LL         = Basis->LL;
   int LR         = Basis->LR;
   int RR         = Basis->RR;
   int RL         = Basis->RL;
   int dim_RR     = Basis->dim_RR;
   int dim_onsite = Basis->dim_onsite;
   int t_elem_num = *elem_num;
   int inv_check,sign,col_LR,col_RR,inv;
   long inv_sup;
   double val,val_LR,val_RR;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   if (strcmp(Sign_Flag, "Yes") == 0) {
      if (strcmp(Type, "LR_RR") == 0) {
         if (Ele_LR[LR]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (strcmp(Type, "RR_LR") == 0) {
         if (Ele_LR[LR]%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         printf("Error in DMRG_MAKE_ELEM_LRRR_LLLRRRRL\n");
         exit(1);
      }
   }
   else {
      sign = 1;
   }
   
   for (i = M_LR->Row[LR]; i < M_LR->Row[LR + 1]; i++) {
      col_LR = M_LR->Col[i];
      val_LR = M_LR->Val[i];
      for (j = M_RR->Row[RR]; j < M_RR->Row[RR + 1]; j++) {
         col_RR    = M_RR->Col[j];
         val_RR    = M_RR->Val[j];
         inv_sup   = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)col_LR*dim_RR*dim_onsite + col_RR*dim_onsite + RL;
         inv       = Inv_LLLRRRRL[inv_sup];
         inv_check = A_Basis->Inv_LLLRRRRL[inv];
         val       = sign*coeef*val_LR*val_RR;
         if (fabs(val) > zero && inv >= 0) {
            if (inv_check == -1) {
               A_Basis->LL_LLLRRRRL[t_elem_num]  = LL;
               A_Basis->LR_LLLRRRRL[t_elem_num]  = col_LR;
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

