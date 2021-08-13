#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_NB_SZBASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   //NB --> Number operator for the Kondo singlet
   //Only for the S = 1/2 KLM
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_NB_SZBASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   if (spin_loc != 1) {
      M->row_dim = 0;
      M->col_dim = 0;
      return;
   }
   
   int dim_charge  = 4;
   int dim_spin    = spin_loc + 1;
   int dim         = dim_spin*dim_charge;
   double zero     = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num;
   double val;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               
               if ((row_c == 1 && col_c == 1) && (row_s == 1 && col_s == 1)) {
                  val = 1.0;
               }
               else if ((row_c == 2 && col_c == 2) && (row_s == 0 && col_s == 0)) {
                  val = 1.0;
               }
               else if ((row_c == 1 && col_c == 2) && (row_s == 1 && col_s == 0)) {
                  val = -1.0;
               }
               else if ((row_c == 2 && col_c == 1) && (row_s == 0 && col_s == 1)) {
                  val = -1.0;
               }
               else {
                  val = 0.0;
               }
               
               val = 0.5*val*coeef;
               
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
