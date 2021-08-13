#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SXL_SZBASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SXL_SZBASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   double zero    = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num,a,b;
   double val,val1,val2;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               a = row_s + 1;
               b = col_s + 1;
               
               if (a == b + 1) {
                  val1 = sqrt((spin_loc*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val1 = 0;
               }
               
               if (a + 1 == b) {
                  val2 = sqrt((spin_loc*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val2 = 0;
               }
               
               val = 0.5*(val1 + val2);
               
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
