#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

void MAKE_SC_MAT_BASIS_1DTKLM_VF(int sc_sz, int ele_1, int ele_2, SC_MAT_1DTKLM_VF *Sc_Mat, MODEL_1DTKLM_VF *Model) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int dim_ccsl  = Model->dim_ccsl_onsite;
   int dim_lspin = Model->dim_lspin;
   int num_ls,num_c,comp_sz,comp_ele_1,comp_ele_2,temp_mat_dim_ccsl,ls_row1,ls_col1,num;
   char C_Name1[100], C_Name2[100];
   
   //Onsite SC correlations
   temp_mat_dim_ccsl = 0;
   for (num = 0; num < dim_ccsl; num++) {
      num_ls  = num/6;
      num_c   = num%6;
      ls_row1 = num_ls/dim_lspin;
      ls_col1 = num_ls%dim_lspin;
      
      if (num_c == 0) {
         //C_1(down)C_1(up)
         comp_ele_1 = 2;
         comp_ele_2 = 0;
         comp_sz    = (ls_col1 - ls_row1)*2;
         sprintf(C_Name1, "Down1");
         sprintf(C_Name2, "Up1");
      }
      else if (num_c == 1) {
         //C_1(down)C_2(down)
         comp_ele_1 = 1;
         comp_ele_2 = 1;
         comp_sz    = (ls_col1 - ls_row1)*2 + 2;
         sprintf(C_Name1, "Down1");
         sprintf(C_Name2, "Down2");
      }
      else if (num_c == 2) {
         //C_1(down)C_2(up)
         comp_ele_1 = 1;
         comp_ele_2 = 1;
         comp_sz    = (ls_col1 - ls_row1)*2;
         sprintf(C_Name1, "Down1");
         sprintf(C_Name2, "Up2");
      }
      else if (num_c == 3) {
         //C_1(up)C_2(down)
         comp_ele_1 = 1;
         comp_ele_2 = 1;
         comp_sz    = (ls_col1 - ls_row1)*2;
         sprintf(C_Name1, "Up1");
         sprintf(C_Name2, "Down2");
      }
      else if (num_c == 4) {
         //C_1(up)C_2(up)
         comp_ele_1 = 1;
         comp_ele_2 = 1;
         comp_sz    = (ls_col1 - ls_row1)*2 - 2;
         sprintf(C_Name1, "Up1");
         sprintf(C_Name2, "Up2");
      }
      else if (num_c == 5) {
         //C_2(down)C_2(up)
         comp_ele_1 = 0;
         comp_ele_2 = 2;
         comp_sz    = (ls_col1 - ls_row1)*2;
         sprintf(C_Name1, "Down2");
         sprintf(C_Name2, "Up2");
      }
      else {
         printf("Error in MAKE_SC_MAT_BASIS_1DTKLM_VF\n");
         exit(1);
      }
      
      
      if (comp_sz == sc_sz && comp_ele_1 == ele_1 && comp_ele_2 == ele_2) {
         Sc_Mat->CCSL_Num[temp_mat_dim_ccsl] = num;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_ccsl], "D[S%d%d%s%s(%d)]", ls_row1 + 1, ls_col1 + 1, C_Name1, C_Name2, Model->tot_site/2);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_ccsl], "S%d%d%s%s(r)", ls_row1 + 1, ls_col1 + 1, C_Name1, C_Name2);
         temp_mat_dim_ccsl++;
      }
      Sc_Mat->CCSL_Sz[num]    = comp_sz;
      Sc_Mat->CCSL_Ele_1[num] = comp_ele_1;
      Sc_Mat->CCSL_Ele_2[num] = comp_ele_2;
   }
   
   Sc_Mat->dim_ccsl = temp_mat_dim_ccsl;
   
   
}
