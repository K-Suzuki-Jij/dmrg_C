#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SP_SZBASIS_HB(int spin, CRS1 *M, double coeef) {
   
   if (spin <= 0) {
      printf("Error in ONSITE_SP_SZBASIS_HB\n");
      printf("2spin=%d\n",spin);
      exit(1);
   }
   
   int dim = spin + 1;
   double zero = pow(10,-15);
   int row,col,elem_num,a,b;
   double val;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row = 0; row < dim; row++) {
      for (col = 0; col < dim; col++) {
         a = row + 1;
         b = col + 1;
         
         if (a + 1 == b) {
            val = sqrt((spin*0.5 + 1)*(a + b - 1.0) - a*b);
         }
         else {
            val = 0;
         }
         
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
