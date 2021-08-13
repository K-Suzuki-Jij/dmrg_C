#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"
#include "onsite.h"

void ONSITE_CCSL_TKLM(int num, int spin_loc, CRS1 *M, double coeef) {
   
   //CCSL(num)
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("ONSITE_CCSL_TKLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge = 4;
   int dim_spin   = spin_loc + 1;
   int dim        = dim_spin*dim_charge*dim_charge;
   int num_ls     = num/6;
   int num_c      = num%6;
   
   CRS1 *M_LS = GET_CRS1(dim, dim*dim);
   CRS1 *M_C1 = GET_CRS1(dim, dim*dim);
   CRS1 *M_C2 = GET_CRS1(dim, dim*dim);
   CRS1 *M_CC = GET_CRS1(dim, dim*dim);
   
   ONSITE_SL_TKLM(num_ls, spin_loc, M_LS, coeef);
   
   if (num_c == 0) {
      //C_1(down)C_1(up)
      ONSITE_CDOWN_1_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CUP_1_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);
   }
   else if (num_c == 1) {
      //C_1(down)C_2(down)
      ONSITE_CDOWN_1_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CDOWN_2_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);

   }
   else if (num_c == 2) {
      //C_1(down)C_2(up)
      ONSITE_CDOWN_1_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CUP_2_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);

   }
   else if (num_c == 3) {
      //C_1(up)C_2(down)
      ONSITE_CUP_1_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CDOWN_2_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);

   }
   else if (num_c == 4) {
      //C_1(up)C_2(up)
      ONSITE_CUP_1_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CUP_2_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);
   }
   else if (num_c == 5) {
      //C_2(down)C_2(up)
      ONSITE_CDOWN_2_SZBASIS_TKLM(spin_loc, M_C1, 1.0);
      ONSITE_CUP_2_SZBASIS_TKLM(spin_loc, M_C2, 1.0);
      MATRIX_PRODUCT_CRS1(M_C1, M_C2, M_CC);
   }
   else {
      printf("Error in ONSITE_CCSL_TKLM\n");
      exit(1);
   }
   
   
   MATRIX_PRODUCT_CRS1(M_CC, M_LS, M);
   FREE_CRS1(M_LS);
   FREE_CRS1(M_C1);
   FREE_CRS1(M_C2);
   FREE_CRS1(M_CC);
   
}

