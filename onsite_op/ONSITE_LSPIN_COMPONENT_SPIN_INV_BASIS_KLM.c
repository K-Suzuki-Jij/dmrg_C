#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_LSPIN_COMPONENT_SPIN_INV_BASIS_KLM(int num, int spin_loc, CRS1 *M, double coeef) {
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_LSPIN_COMPONENT_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   int sc_row     = num/dim_spin;
   int sc_col     = num%dim_spin;
   double zero    = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num;
   double val;

   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
             
               if (row_s == sc_row && col_s == sc_col) {
                  val = coeef;
               }
               else {
                  val = 0.0;
               }
               
               if (fabs(val) > zero && row_c == col_c) {
                  M->Val[elem_num] = val;
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

