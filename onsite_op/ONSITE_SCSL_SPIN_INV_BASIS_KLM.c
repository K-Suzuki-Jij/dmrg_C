#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void ONSITE_SCSL_SPIN_INV_BASIS_KLM(int spin_loc, CRS1 *M, double coeef) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   
   if (spin_loc <= 0) {
      printf("Error in ONSITE_SCSL_SPIN_INV_BASIS_KLM\n");
      printf("2spin=%d\n",spin_loc);
      exit(1);
   }
   
   int dim_charge  = 4;
   int dim_spin    = spin_loc + 1;
   int dim         = dim_spin*dim_charge;
   int dim_parity  = 2;
   double zero     = pow(10,-15);
   int row_c,col_c,row_s,col_s,elem_num,mz_row_s,mz_col_s,parity_row_s,parity_col_s,parity_row_c,parity_col_c;
   double val,val_szc,val_spc,val_smc,val_szl,val_spl,val_sml,c_p,c_m;
   
   M->Row[0] = 0;
   elem_num = 0;
   for (row_c = 0; row_c < dim_charge; row_c++) {
      for (row_s = 0; row_s < dim_spin; row_s++) {
         for (col_c = 0; col_c < dim_charge; col_c++) {
            for (col_s = 0; col_s < dim_spin; col_s++) {
               parity_row_c = row_c/dim_parity;
               parity_col_c = col_c/dim_parity;
               
               //SzC
               if (parity_row_c != parity_col_c) {
                  val_szc = 0.5;
               }
               else {
                  val_szc = 0.0;
               }
               
               //SpC
               if (col_c == 1 && row_c == 1) {
                  val_spc = 0.5;
               }
               else if (col_c == 2 && row_c == 2) {
                  val_spc = -0.5;
               }
               else if (col_c == 2 && row_c == 1) {
                  val_spc = -0.5;
               }
               else if (col_c == 1 && row_c == 2) {
                  val_spc = 0.5;
               }
               else {
                  val_spc = 0.0;
               }
               
               //SmC
               if (col_c == 1 && row_c == 1) {
                  val_smc = 0.5;
               }
               else if (col_c == 2 && row_c == 2) {
                  val_smc = -0.5;
               }
               else if (col_c == 2 && row_c == 1) {
                  val_smc = 0.5;
               }
               else if (col_c == 1 && row_c == 2) {
                  val_smc = -0.5;
               }
               else {
                  val_smc = 0.0;
               }
            
               if (row_c == 0 || row_c == 3 || col_c == 0 || col_c == 3) {
                  val_szc = 0.0;
                  val_spc = 0.0;
                  val_smc = 0.0;
               }
               
               ///////SL
               mz_row_s     = spin_loc - 2*(row_s/dim_parity);
               parity_row_s = row_s%dim_parity;
               mz_col_s     = spin_loc - 2*(col_s/dim_parity);
               parity_col_s = col_s%dim_parity;
               
               //SzL
               if (mz_row_s == mz_col_s && parity_row_s != parity_col_s) {
                  val_szl = 0.5*mz_col_s;
               }
               else {
                  val_szl = 0;
               }
               
               //SpL
               c_p = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col_s*0.5*(mz_col_s*0.5 + 1.0)));
               c_m = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col_s*0.5*(mz_col_s*0.5 - 1.0)));
               
               if (mz_col_s > 2) {
                  if (mz_col_s + 2 == mz_row_s) {
                     val_spl = 0.5*c_p;
                  }
                  else if (mz_col_s - 2 == mz_row_s) {
                     if (parity_col_s == parity_row_s) {
                        val_spl = 0.5*c_m;
                     }
                     else {
                        val_spl = -0.5*c_m;
                     }
                  }
                  else {
                     val_spl = 0.0;
                  }
               }
               else if (mz_col_s == 2) {
                  if (mz_col_s + 2 == mz_row_s) {
                     val_spl = 0.5*c_p;
                  }
                  else if (mz_col_s - 2 == mz_row_s) {
                     if (parity_col_s == 0) {
                        val_spl = c_m/sqrt(2);
                     }
                     else {
                        val_spl = -c_m/sqrt(2);
                     }
                  }
                  else {
                     val_spl = 0.0;
                  }
               }
               else if (mz_col_s == 1) {
                  if (mz_col_s + 2 == mz_row_s) {
                     val_spl = 0.5*c_p;
                  }
                  else if (2 - mz_col_s == mz_row_s) {
                     if (parity_col_s == 0) {
                        val_spl = 0.5*c_m;
                     }
                     else {
                        val_spl = -0.5*c_m;
                     }
                  }
                  else {
                     val_spl = 0.0;
                  }
               }
               else if (mz_col_s == 0) {
                  if (mz_col_s + 2 == mz_row_s) {
                     val_spl = c_p/sqrt(2);
                  }
                  else {
                     val_spl = 0.0;
                  }
               }
               else {
                  printf("Error in ONSITE_SCSL_SPIN_INV_BASIS_KLM\n");
                  exit(1);
               }
               
               //SmL
               c_p = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col_s*0.5*(mz_col_s*0.5 - 1.0)));
               c_m = sqrt(fabs(spin_loc*0.5*(spin_loc*0.5 + 1.0) - mz_col_s*0.5*(mz_col_s*0.5 + 1.0)));
               
               if (mz_col_s > 2) {
                  if (mz_col_s - 2 == mz_row_s) {
                     val_sml = 0.5*c_p;
                  }
                  else if (mz_col_s + 2 == mz_row_s) {
                     if (parity_col_s == parity_row_s) {
                        val_sml = 0.5*c_m;
                     }
                     else {
                        val_sml = -0.5*c_m;
                     }
                  }
                  else {
                     val_sml = 0.0;
                  }
               }
               else if (mz_col_s == 2) {
                  if (2 - mz_col_s == mz_row_s) {
                     val_sml = c_p/sqrt(2);
                  }
                  else if (mz_col_s + 2 == mz_row_s) {
                     if (parity_col_s == parity_row_s) {
                        val_sml = 0.5*c_m;
                     }
                     else {
                        val_sml = -0.5*c_m;
                     }
                  }
                  else {
                     val_sml = 0.0;
                  }
               }
               else if (mz_col_s == 1) {
                  if (2 - mz_col_s == mz_row_s) {
                     if (parity_row_s == 0) {
                        val_sml = 0.5*c_p;
                     }
                     else {
                        val_sml = -0.5*c_p;
                     }
                  }
                  else if (mz_col_s + 2 == mz_row_s) {
                     if (parity_col_s == parity_row_s) {
                        val_sml = 0.5*c_m;
                     }
                     else {
                        val_sml = -0.5*c_m;
                     }
                  }
                  else {
                     val_sml = 0.0;
                  }
               }
               else if (mz_col_s == 0) {
                  if (2 - mz_col_s == mz_row_s) {
                     if (parity_row_s == 0) {
                        val_sml = c_p/sqrt(2);
                     }
                     else {
                        val_sml = -c_p/sqrt(2);
                     }
                  }
                  else {
                     val_sml = 0.0;
                  }
               }
               else {
                  printf("Error in ONSITE_SCSL_SPIN_INV_BASIS_KLM\n");
                  exit(1);
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
