#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

void MAKE_LSC_MAT_BASIS_1DKLM_VF(int sc_sz, SC_MAT_1DKLM_VF *Sc_Mat, int lspin, int tot_site) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////

   int dim_lspin = lspin + 1;
   int dim_lscc  = dim_lspin*dim_lspin;
   int dim_lsc   = dim_lspin*dim_lspin*2;
   int num1,num2,ls_row1,ls_row2,ls_col1,ls_col2,comp_sz1,comp_sz2,temp_mat_dim_lscc_1,temp_mat_dim_lscc_2,temp_mat_dim_lsc_lsc,c_op1,c_op2;
   char C_Name1[100], C_Name2[100];
   
   //Onsite SC correlations S*Up(i)*Down(i)
   temp_mat_dim_lscc_1 = 0;
   for (num1 = 0; num1 < dim_lscc; num1++) {
      ls_row1 = num1/dim_lspin;
      ls_col1 = num1%dim_lspin;
      
      comp_sz1 = (ls_col1 - ls_row1)*2;
      
      if (comp_sz1 == sc_sz) {
         Sc_Mat->CC_1_Num[temp_mat_dim_lscc_1] = num1;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_lscc_1], "D[S%d%dDownUp(%d)]", ls_row1 + 1, ls_col1 + 1, tot_site/2);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_lscc_1], "S%d%dDownUp(r)", ls_row1 + 1, ls_col1 + 1);
         temp_mat_dim_lscc_1++;
      }
      Sc_Mat->CC_Sz[num1] = comp_sz1;
   }
   
   //Intersite SC correlations S1*C_1(i)S2*C_2(i+1)
   temp_mat_dim_lsc_lsc = 0;
   for (num1 = 0; num1 < dim_lsc; num1++) {
      c_op1   = num1%2 + 1;
      ls_row1 = (num1/2)/dim_lspin;
      ls_col1 = (num1/2)%dim_lspin;
      if (c_op1 == 1) {
         comp_sz1 = (ls_col1 - ls_row1)*2 - 1;
         sprintf(C_Name1, "Up");
      }
      else if (c_op1 == 2) {
         comp_sz1 = (ls_col1 - ls_row1)*2 + 1;
         sprintf(C_Name1, "Down");
      }
      else {
         printf("Error in MAKE_LSC_MAT_BASIS_1DKLM_VF\n");
         exit(1);
      }
      for (num2 = 0; num2 < dim_lsc; num2++) {
         c_op2   = num2%2 + 1;
         ls_row2 = (num2/2)/dim_lspin;
         ls_col2 = (num2/2)%dim_lspin;
         if (c_op2 == 1) {
            comp_sz2 = (ls_col2 - ls_row2)*2 - 1;
            sprintf(C_Name2, "Up");
         }
         else if (c_op2 == 2) {
            comp_sz2 = (ls_col2 - ls_row2)*2 + 1;
            sprintf(C_Name2, "Down");
         }
         else {
            printf("Error in MAKE_LSC_MAT_BASIS_1DKLM_VF\n");
            exit(1);
         }
         if (comp_sz1 + comp_sz2 == sc_sz) {
            Sc_Mat->C_Num1[temp_mat_dim_lsc_lsc] = num1;
            Sc_Mat->C_Num2[temp_mat_dim_lsc_lsc] = num2;
            sprintf(Sc_Mat->Row_Name[temp_mat_dim_lsc_lsc + temp_mat_dim_lscc_1], "D[S%d%d%s(%d)S%d%d%s(%d)]", ls_row1 + 1, ls_col1 + 1, C_Name1, tot_site/2, ls_row2 + 1, ls_col2 + 1, C_Name2, tot_site/2 - 1);
            sprintf(Sc_Mat->Col_Name[temp_mat_dim_lsc_lsc + temp_mat_dim_lscc_1], "S%d%d%s(r)S%d%d%s(r+1)", ls_row1 + 1, ls_col1 + 1, C_Name1, ls_row2 + 1, ls_col2 + 1, C_Name2);
            temp_mat_dim_lsc_lsc++;
         }
      }
      Sc_Mat->C_Sz[num1] = comp_sz1;
   }
   
   temp_mat_dim_lscc_2 = 0;
   //Onsite SC correlations S*Even*Odd(i+1)
   for (num1 = 0; num1 < dim_lscc; num1++) {
      ls_row1 = num1/dim_lspin;
      ls_col1 = num1%dim_lspin;
      comp_sz1 = (ls_col1 - ls_row1)*2;
      if (comp_sz1 == sc_sz) {
         Sc_Mat->CC_2_Num[temp_mat_dim_lscc_2] = num1;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_lsc_lsc + temp_mat_dim_lscc_1 + temp_mat_dim_lscc_2], "D[S%d%dDownUp(%d)]", ls_row1 + 1, ls_col1 + 1, tot_site/2 - 1);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_lsc_lsc + temp_mat_dim_lscc_1 + temp_mat_dim_lscc_2], "S%d%dDownUp(r+1)", ls_row1 + 1, ls_col1 + 1);
         temp_mat_dim_lscc_2++;
      }
   }
   
   Sc_Mat->dim_cc_1 = temp_mat_dim_lscc_1;
   Sc_Mat->dim_cc_2 = temp_mat_dim_lscc_2;
   Sc_Mat->dim_c_c  = temp_mat_dim_lsc_lsc;
   Sc_Mat->dim_tot  = temp_mat_dim_lscc_1 + temp_mat_dim_lscc_2 + temp_mat_dim_lsc_lsc;
   
   //Check Point
   if (temp_mat_dim_lscc_1 != temp_mat_dim_lscc_2) {
      printf("Error in MAKE_LSC_MAT_BASIS_1DKLM_VF\n");
      exit(1);
   }
   
}
