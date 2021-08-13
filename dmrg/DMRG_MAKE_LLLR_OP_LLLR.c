#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dmrg.h"

void DMRG_MAKE_LLLR_OP_LLLR(CRS1 *M_LL, CRS1 *M_LR, double coeef, int *Ele_LL, char Sign[], char Type[], DMRG_BASIS_LLLR *Dmrg_Basis, CRS1 *Out) {
   
   int dim_LLLR = Dmrg_Basis->dim_LLLR;
   long i,j,basis,elem_num;
   int LL,LR,inv,sign,A_LL,A_LR;
   double val_LL,val_LR;
   
   elem_num = 0;
   for (basis = 0; basis < dim_LLLR; basis++) {
      LL = Dmrg_Basis->LL_LLLR[basis];
      LR = Dmrg_Basis->LR_LLLR[basis];
      
      if (strcmp(Sign, "Yes") == 0) {
         if (strcmp(Type, "LL_LR") == 0) {
            if (Ele_LL[LL]%2 == 0) {
               sign = -1;
            }
            else {
               sign = 1;
            }
         }
         else if (strcmp(Type, "LR_LL") == 0) {
            if (Ele_LL[LL]%2 == 0) {
               sign = 1;
            }
            else {
               sign = -1;
            }
         }
         else {
            printf("Error in DMRG_MAKE_LLLR_OP_LLLR\n");
            exit(1);
         }
      }
      else {
         sign = 1;
      }
      
      for (i = M_LL->Row[LL]; i < M_LL->Row[LL + 1]; i++) {
         A_LL   = M_LL->Col[i];
         val_LL = M_LL->Val[i];
         for (j = M_LR->Row[LR]; j < M_LR->Row[LR + 1]; j++) {
            A_LR   = M_LR->Col[j];
            val_LR = M_LR->Val[j];
            inv = Dmrg_Basis->Inv_LLLR[A_LL][A_LR];
            if (inv >= 0) {
               //Check Point
               if (elem_num >= Out->max_val) {
                  printf("Error in DMRG_MAKE_LLLR_OP_LLLR\n");
                  printf("Need more Out->max_val = %ld\n", Out->max_val);
                  exit(1);
               }
               Out->Val[elem_num] = sign*coeef*val_LL*val_LR;
               Out->Col[elem_num] = inv;
               elem_num++;
            }
         }
      }
      
      //Check Point
      if (basis + 1 >= Out->max_row) {
         printf("Error in DMRG_MAKE_LLLR_OP_LLLR\n");
         printf("Need more Out->max_row = %d\n", Out->max_row);
         exit(1);
      }
      
      Out->Row[basis + 1] = elem_num;
      
   }
   
   Out->row_dim = dim_LLLR;
   Out->col_dim = dim_LLLR;
   
   SORT_COLUMN_CRS1(Out,1);
   
}
