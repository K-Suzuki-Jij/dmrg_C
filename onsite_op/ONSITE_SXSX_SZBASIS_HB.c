#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SXSX_SZBASIS_HB(int spin, CRS1 *M, double coeef) {
   
   if (spin <= 0) {
      printf("Error in ONSITE_SXSX_SZBASIS_HB\n");
      printf("2spin=%d\n",spin);
      exit(1);
   }
   
   int dim = spin + 1;
   double zero = pow(10,-15);
   int row,col,elem_num,a,b;
   double val,val1,val2,val3,val4;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row = 0; row < dim; row++) {
      for (col = 0; col < dim; col++) {
         a = row + 1;
         b = col + 1;
         
         if (a == b + 2) {
            val1 = sqrt((spin*0.5 + 1.0)*(2.0*a - 2.0) - a*(a - 1.0))*sqrt((spin*0.5 + 1.0)*(2.0*a - 4.0) - (a - 1.0)*(a - 2.0));
         }
         else {
            val1 = 0;
         }
         
         if (a == b) {
            val2 = (spin*0.5 + 1.0)*(2.0*a - 2.0) - a*(a - 1.0);
            val3 = (spin*0.5 + 1.0)*(2.0*a) - a*(a + 1.0);
         }
         else {
            val2 = 0;
            val3 = 0;
         }
         
         if (a + 2 == b) {
            val4 = sqrt((spin*0.5 + 1.0)*(2.0*a) - a*(a + 1.0))*sqrt((spin*0.5 + 1.0)*(2.0*a + 2.0) - (a + 2.0)*(a + 1.0));
         }
         else {
            val4 = 0;
         }

         val  = 0.25*(val1 + val2 + val3 + val4);
         
         val = val*coeef;
         
         if (fabs(val) > zero) {
            M->Val[elem_num] = val;
            M->Col[elem_num] = col;
            elem_num++;
         }
      }
      M->Row[row + 1] = elem_num;
   }
   
   M->row_dim = dim;
   M->col_dim = dim;
   
}

