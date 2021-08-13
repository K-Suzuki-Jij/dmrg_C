#include <stdio.h>
#include <stdlib.h>
#include "SML.h"

void ONSITE_DIAG_SPIN_INV_BASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   if (spin_loc <= 0) {
      printf("ONSITE_DIAG_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   int row_c,col_c,row_s,col_s,elem_num;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               
               if (row_c == col_c && row_s == col_s) {
                  M->Val[elem_num] = coeef;
                  M->Col[elem_num] = col_c*dim_spin + col_s;
                  elem_num++;
               }
            }
         }
         M->Row[row_c*dim_spin + row_s + 1] = elem_num;
      }
   }
   
   M->row_dim = dim;
   M->col_dim = dim;
   
}
