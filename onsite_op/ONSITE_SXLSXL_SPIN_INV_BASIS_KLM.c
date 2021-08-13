#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"
#include "onsite.h"

void ONSITE_SXLSXL_SPIN_INV_BASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SXLSXL_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   
   CRS1 *Sx = GET_CRS1(dim, dim*dim);
   ONSITE_SXL_SPIN_INV_BASIS_KLM(spin_loc, Sx, 1.0);
   MATRIX_PRODUCT_CRS1(Sx, Sx, M);
   
   int row;
   long iter;
   
   for (row = 0; row < M->row_dim; row++) {
      for (iter = M->Row[row]; iter < M->Row[row + 1]; iter++) {
         M->Val[iter] = M->Val[iter]*coeef;
      }
   }
   
   FREE_CRS1(Sx);
   
}
