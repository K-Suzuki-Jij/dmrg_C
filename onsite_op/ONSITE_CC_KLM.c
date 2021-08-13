#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_CC_KLM(int num, int spin_loc, CRS1 *M, double coeef) {
      
   //CC(num)
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0 || num != 0) {
      printf("ONSITE_CC_KLM\n");
      printf("2spin=%d,num=%d\n",spin_loc,num);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   double zero    = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num;
   double val;
   
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               
               if (row_s == col_s && col_c == 3 && row_c == 0) {
                  val = coeef;
               }
               else {
                  val = 0.0;
               }
               
               if (fabs(val) > zero) {
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

