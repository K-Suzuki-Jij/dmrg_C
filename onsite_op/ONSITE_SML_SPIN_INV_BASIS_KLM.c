#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"
#include "onsite.h"

void ONSITE_SML_SPIN_INV_BASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SPL_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge;
   
   CRS1 *Sp = GET_CRS1(dim, dim*dim);
   ONSITE_SPL_SPIN_INV_BASIS_KLM(spin_loc, Sp, coeef);
   MATRIX_TRANSPOSE_CRS1(Sp, M);
   FREE_CRS1(Sp);
    
   
}
