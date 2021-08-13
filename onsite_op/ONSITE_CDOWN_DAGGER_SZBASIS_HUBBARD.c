#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_CDOWN_DAGGER_SZBASIS_HUBBARD(CRS1 *M, double coeef) {
      
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int dim     = 4;
   double zero = pow(10,-15);
   int row,col,elem_num;
   double val;
   
   M->Row[0] = 0;
   elem_num = 0;
   
   for (row = 0; row < dim; row++) {
      for (col = 0; col < dim; col++) {
         
         if (col == 0 && row == 2) {
            val = 1.0;
         }
         else if (col == 1 && row == 3) {
            val = -1.0;
         }
         else {
            val = 0.0;
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
