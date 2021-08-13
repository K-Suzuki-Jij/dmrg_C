#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LRRL_LRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LRRL, CRS1 *M_LR, CRS1 *M_RL, double coeef, int *Ele_LR, int *elem_num, char Type[], char Sign_Flag[]) {
   
   long i,j;
   
   int LR = Basis->LR;
   int RL = Basis->RL;
   int t_elem_num = *elem_num;
   int inv_check,sign,col_LR,col_RL,inv;
   double val,val_LR,val_RL;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   if (strcmp(Sign_Flag, "Yes") == 0) {
      if (strcmp(Type, "LR_RL") == 0) {
         if (Ele_LR[LR]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (strcmp(Type, "RL_LR") == 0) {
         if (Ele_LR[LR]%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         printf("Error in DMRG_MAKE_ELEM_LRRL_LRRL\n");
         exit(1);
      }
   }
   else {
      sign = 1;
   }
   
   for (i = M_LR->Row[LR]; i < M_LR->Row[LR + 1]; i++) {
      col_LR = M_LR->Col[i];
      val_LR = M_LR->Val[i];
      for (j = M_RL->Row[RL]; j < M_RL->Row[RL + 1]; j++) {
         col_RL    = M_RL->Col[j];
         val_RL    = M_RL->Val[j];
         inv       = Inv_LRRL[col_LR][col_RL];
         inv_check = A_Basis->Inv_LLLRRRRL[inv];
         val       = sign*coeef*val_LR*val_RL;
         if (fabs(val) > zero && inv >= 0) {
            if (inv_check == -1) {
               A_Basis->LR_LLLRRRRL[t_elem_num]  = col_LR;
               A_Basis->RL_LLLRRRRL[t_elem_num]  = col_RL;
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
