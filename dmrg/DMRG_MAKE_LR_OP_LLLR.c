#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dmrg.h"

void DMRG_MAKE_LR_OP_LLLR(CRS1 *M_On, int *Ele_LL, char Sign[], DMRG_BASIS_LLLR *Dmrg_Basis, CRS1 *Out) {
   
   int dim_LLLR   = Dmrg_Basis->dim_LLLR;
   long iter,basis,elem_num;
   int LL,LR,inv,sign;
   
   elem_num = 0;
   for (basis = 0; basis < dim_LLLR; basis++) {
      LL = Dmrg_Basis->LL_LLLR[basis];
      LR = Dmrg_Basis->LR_LLLR[basis];
      
      if (strcmp(Sign, "Yes") == 0) {
         if (Ele_LL[LL]%2 == 0) {
            sign = -1;
         }
         else {
            sign = 1;
         }
      }
      else {
         sign = 1;
      }
      
      for (iter = M_On->Row[LR]; iter < M_On->Row[LR + 1]; iter++) {
         inv = Dmrg_Basis->Inv_LLLR[LL][M_On->Col[iter]];
         if (inv >= 0) {
            //Check Point
            if (elem_num >= Out->max_val) {
               printf("Error in DMRG_MAKE_LR_OP_LLLR\n");
               printf("Need more Out->max_val = %ld\n", Out->max_val);
               exit(1);
            }
            Out->Val[elem_num] = sign*M_On->Val[iter];
            Out->Col[elem_num] = inv;
            elem_num++;
         }
      }
      
      //Check Point
      if (basis + 1 >= Out->max_row) {
         printf("Error in DMRG_MAKE_LR_OP_LLLR\n");
         printf("Need more Out->max_row = %d\n", Out->max_row);
         exit(1);
      }
      
      Out->Row[basis + 1] = elem_num;
      
   }
   
   Out->row_dim = dim_LLLR;
   Out->col_dim = dim_LLLR;
   
   SORT_COLUMN_CRS1(Out,1);
   
}
