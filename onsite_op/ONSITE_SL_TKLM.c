#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SL_TKLM(int num, int spin_loc, CRS1 *M, double coeef) {
   
   //CCSL(num)
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("ONSITE_SL_TKLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge*dim_charge;
   int ls_row     = num/dim_spin;
   int ls_col     = num%dim_spin;
   double zero    = pow(10,-15);
   int row_c_1,col_c_1,row_c_2,col_c_2,row_s,col_s,elem_num;
   double val;
   
   
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c_1 = 0; row_c_1 < dim_charge; row_c_1++) {
      for (row_c_2 = 0; row_c_2 < dim_charge; row_c_2++) {
         for (row_s = 0; row_s < dim_spin; row_s++) {
            for (col_c_1 = 0; col_c_1 < dim_charge; col_c_1++) {
               for (col_c_2 = 0; col_c_2 < dim_charge; col_c_2++) {
                  for (col_s = 0; col_s < dim_spin; col_s++) {
                     
                     if (row_s == ls_row && col_s == ls_col && col_c_1 == row_c_1 && row_c_2 == col_c_2) {
                        val = coeef;
                     }
                     else {
                        val = 0.0;
                     }
                     
                     if (fabs(val) > zero) {
                        M->Val[elem_num] = val;
                        M->Col[elem_num] = col_c_1*dim_charge*dim_spin + col_c_2*dim_spin + col_s;
                        elem_num++;
                     }
                  }
               }
            }
            M->Row[row_c_1*dim_charge*dim_spin + row_c_2*dim_spin + row_s + 1] = elem_num;
         }
      }
   }
   
   M->row_dim = dim;
   M->col_dim = dim;
   
}

