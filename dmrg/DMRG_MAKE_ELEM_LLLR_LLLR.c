#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmrg.h"

void DMRG_MAKE_ELEM_LLLR_LLLR(DMRG_BASIS_ONSITE *Basis, DMRG_A_BASIS *A_Basis, int **Inv_LLLR, CRS1 *M_LL, CRS1 *M_LR, double coeef, int *LL_Ele, int *elem_num, char Type[], char Sign_Flag[]) {
   
   long i,j;
   
   int LL = Basis->LL;
   int LR = Basis->LR;
   
   int t_elem_num = *elem_num;
   int inv_check,sign,col_LL,col_LR,inv;
   
   double val,val_LL,val_LR;
   double zero = 0.000000000000001;//pow(10,-15)
   
   if (fabs(coeef) < zero) {
      return;
   }
   
   if (strcmp(Sign_Flag, "Yes") == 0) {
      if (strcmp(Type, "LL_LR") == 0) {
         if (LL_Ele[LL]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else if (strcmp(Type, "LR_LL") == 0) {
         if (LL_Ele[LL]%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         printf("Error in DMRG_MAKE_ELEM_LLLR_LLLR\n");
         exit(1);
      }
   }
   else {
      sign = 1;
   }
   
   for (i = M_LL->Row[LL]; i < M_LL->Row[LL + 1]; i++) {
      col_LL = M_LL->Col[i];
      val_LL = M_LL->Val[i];
      for (j = M_LR->Row[LR]; j < M_LR->Row[LR + 1]; j++) {
         col_LR    = M_LR->Col[j];
         val_LR    = M_LR->Val[j];
         inv       = Inv_LLLR[col_LL][col_LR];
         inv_check = A_Basis->Inv_LLLRRRRL[inv];
         val       = sign*coeef*val_LL*val_LR;
         if (fabs(val) > zero && inv >= 0) {
            if (inv_check == -1) {
               A_Basis->LL_LLLRRRRL[t_elem_num]  = col_LL;
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
   }
   
   *elem_num = t_elem_num;
   
}
