#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_CDOWN_2_SZBASIS_TKLM(int spin_loc, CRS1 *M, double coeef) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_CDOWN_2_SZBASIS_TKLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge  = 4;
   int dim_spin    = spin_loc + 1;
   int dim         = dim_spin*dim_charge*dim_charge;
   double zero     = pow(10,-15);
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
                     
                     if (col_c_2 == 2 && row_c_2 == 0) {
                        val = 1.0;
                     }
                     else if (col_c_2 == 3 && row_c_2 == 1) {
                        val = -1.0;
                     }
                     else {
                        val = 0.0;
                     }
                     
                     if (row_c_1 == 1 || row_c_1 == 2) {
                        val = -1.0*val;
                     }
                     
                     val = val*coeef;
                     
                     if (fabs(val) > zero && row_s == col_s && row_c_1 == col_c_1) {
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
