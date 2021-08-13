#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SXL_SPIN_INV_BASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   //|S,mz,ev> |S,mz,od> |S,mz-1,ev> |S,mz-1,od> ...
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SXL_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   int dim_parity = 2;
   double zero    = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num,mz_row,mz_col,parity_row,parity_col;
   double val,c_p,c_m;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               mz_row     = spin_loc - 2*(row_s/dim_parity);
               parity_row = row_s%dim_parity;
               mz_col     = spin_loc - 2*(col_s/dim_parity);
               parity_col = col_s%dim_parity;
               
               c_p = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col*0.5*(+1.0*mz_col*0.5 + 1.0)));
               c_m = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col*0.5*(+1.0*mz_col*0.5 - 1.0)));
               
               if (parity_col == parity_row) {
                  if (mz_col > 2) {
                     if (mz_col + 2 == mz_row) {
                        val = 0.5*c_p;
                     }
                     else if (mz_col - 2 == mz_row) {
                        val = 0.5*c_m;
                     }
                     else {
                        val = 0.0;
                     }
                  }
                  else if (mz_col == 2) {
                     if (mz_col + 2 == mz_row) {
                        val = 0.5*c_p;
                     }
                     else if (mz_col - 2 == mz_row) {
                        if (parity_col == 0) {
                           val = 0.5*c_m*sqrt(2);
                        }
                        else {
                           val = 0.0;
                        }
                     }
                     else {
                        val = 0.0;
                     }
                  }
                  else if (mz_col == 1) {
                     if (mz_col + 2 == mz_row) {
                        val = 0.5*c_p;
                     }
                     else if (2 - mz_col == mz_row) {
                        if (parity_col == 0) {
                           val = 0.5*c_m;
                        }
                        else {
                           val = -0.5*c_m;
                        }
                     }
                     else {
                        val = 0.0;
                     }
                  }
                  else if (mz_col == 0) {
                     if (mz_col + 2 == mz_row) {
                        if (parity_col == 0) {
                           val = 0.5*c_p*sqrt(2);
                        }
                        else {
                           printf("Error in ONSITE_SXL_SPIN_INV_BASIS_KLM\n");
                           exit(1);
                        }
                     }
                     else {
                        val = 0.0;
                     }
                  }
                  else {
                     printf("Error in ONSITE_SXL_SPIN_INV_BASIS_KLM\n");
                     exit(1);
                  }
               }
               else {
                  val = 0.0;
               }
               
               val = val*coeef;
               
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
