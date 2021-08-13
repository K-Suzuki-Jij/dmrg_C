#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SCSL_SZBASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SCSL_SZBASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge  = 4;
   int dim_spin    = spin_loc + 1;
   int dim         = dim_spin*dim_charge;
   int spin_charge = 1;
   double zero     = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num,a,b;
   double val,val_szc,val_spc,val_smc,val_szl,val_spl,val_sml;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               
               ///////SC
               a = row_c;
               b = col_c;
               
               //SzC
               val_szc  = (spin_charge*0.5 + 1.0 - b)*DELTA_FUNCTION(a, b);
               
               //SpC
               if (a + 1 == b) {
                  val_spc = sqrt((spin_charge*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val_spc = 0.0;
               }
               
               //SmC
               if (a == b + 1) {
                  val_smc = sqrt((spin_charge*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val_smc = 0;
               }
               
               
               if (row_c == 0 || row_c == 3 || col_c == 0 || col_c == 3) {
                  val_szc = 0.0;
                  val_spc = 0.0;
                  val_smc = 0.0;
               }
               
               
               ///////SL
               a = row_s + 1;
               b = col_s + 1;
               
               //SzL
               val_szl  = (spin_loc*0.5 + 1.0 - b)*DELTA_FUNCTION(a, b);
               
               //SpL
               if (a + 1 == b) {
                  val_spl = sqrt((spin_loc*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val_spl = 0;
               }
               
               //SmL
               if (a == b + 1) {
                  val_sml = sqrt((spin_loc*0.5 + 1)*(a + b - 1.0) - a*b);
               }
               else {
                  val_sml = 0;
               }
               
               val = val_szc*val_szl + 0.5*val_spc*val_sml + 0.5*val_smc*val_spl;
               
               val = val*coeef;
               
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
