#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_RRRL_RRRL(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_RRRL, CRS1 *M_RR, CRS1 *M_RL, double coeef, int *Ele_RR, int *elem_num, char Type[], char Sign_Flag[]) {
   
   long i,j;
   
   int RR  = Basis->RR;
   int RL  = Basis->RL;
   int t_elem_num = *elem_num;
   int inv_check,sign,col_RL,col_RR,inv;
   double val,val_RL,val_RR;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   if (strcmp(Sign_Flag, "Yes") == 0) {
      if (strcmp(Type, "RR_RL") == 0) {
         if (Ele_RR[RR]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (strcmp(Type, "RL_RR") == 0) {
         if (Ele_RR[RR]%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         printf("Error in DMRG_MAKE_ELEM_RRRL_RRRL\n");
         exit(1);
      }
   }
   else {
      sign = 1;
   }
   
   for (i = M_RR->Row[RR]; i < M_RR->Row[RR + 1]; i++) {
      col_RR = M_RR->Col[i];
      val_RR = M_RR->Val[i];
      for (j = M_RL->Row[RL]; j < M_RL->Row[RL + 1]; j++) {
         col_RL    = M_RL->Col[j];
         val_RL    = M_RL->Val[j];
         inv       = Inv_RRRL[col_RR][col_RL];
         inv_check = A_Basis->Inv_LLLRRRRL[inv];
         val       = sign*coeef*val_RR*val_RL;
         if (fabs(val) > zero && inv >= 0) {
            if (inv_check == -1) {
               A_Basis->RR_LLLRRRRL[t_elem_num]  = col_RR;
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

