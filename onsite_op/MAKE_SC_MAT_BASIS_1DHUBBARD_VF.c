#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

void MAKE_SC_MAT_BASIS_1DHUBBARD_VF(int sc_sz, SC_MAT_1DHUBBARD_VF *Sc_Mat, int tot_site) {
   
   ///////////////////////////////////////
   // # <->  [Cherge  ] -- (N,  2*sz)
   // 0 <->  [        ] -- (0,  0   )
   // 1 <->  [up      ] -- (1,  1   )
   // 2 <->  [down    ] -- (1, -1   )
   // 3 <->  [up&down ] -- (2,  0   )
   ///////////////////////////////////////
   
   int dim_cc_onsite = 1;
   int dim_c_onsite  = 2;
   int num1,num2,comp_sz1,comp_sz2,temp_mat_dim_cc_1,temp_mat_dim_cc_2,temp_mat_dim_c_c;
   char C_Name1[100], C_Name2[100];
   
   //Onsite SC correlations Down(i)*Up(i)
   temp_mat_dim_cc_1 = 0;
   for (num1 = 0; num1 < dim_cc_onsite; num1++) {
      comp_sz1 = 0;
      if (comp_sz1 == sc_sz) {
         Sc_Mat->CC_1_Num[temp_mat_dim_cc_1] = num1;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_cc_1], "D[DownUp(%d)]", tot_site/2);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_cc_1], "DownUp(r)");
         temp_mat_dim_cc_1++;
      }
      Sc_Mat->CC_Sz[num1] = comp_sz1;
   }
   
   //Intersite SC correlations C_1(i)C_2(i+1)
   temp_mat_dim_c_c = 0;
   for (num1 = 0; num1 < dim_c_onsite; num1++) {
     
      if (num1 == 0) {
         comp_sz1 = -1;
         sprintf(C_Name1, "Up");
      }
      else if (num1 == 1) {
         comp_sz1 = +1;
         sprintf(C_Name1, "Down");
      }
      else {
         printf("Error in MAKE_SC_MAT_BASIS_1DHUBBARD_VF\n");
         exit(1);
      }
      
      for (num2 = 0; num2 < dim_c_onsite; num2++) {
      
         if (num2 == 0) {
            comp_sz2 = -1;
            sprintf(C_Name2, "Up");
         }
         else if (num2 == 1) {
            comp_sz2 = +1;
            sprintf(C_Name2, "Down");
         }
         else {
            printf("Error in MAKE_SC_MAT_BASIS_1DHUBBARD_VF\n");
            exit(1);
         }
         if (comp_sz1 + comp_sz2 == sc_sz) {
            Sc_Mat->C_Num1[temp_mat_dim_c_c] = num1;
            Sc_Mat->C_Num2[temp_mat_dim_c_c] = num2;
            sprintf(Sc_Mat->Row_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1], "D[%s(%d)%s(%d)]", C_Name1, tot_site/2, C_Name2, tot_site/2 - 1);
            sprintf(Sc_Mat->Col_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1], "%s(r)%s(r+1)", C_Name1, C_Name2);
            temp_mat_dim_c_c++;
         }
      }
      Sc_Mat->C_Sz[num1] = comp_sz1;
   }
   
   //Onsite SC correlations Down(i)*Up(i)
   temp_mat_dim_cc_2 = 0;
   for (num2 = 0; num2 < dim_cc_onsite; num2++) {
      comp_sz2 = 0;
      if (comp_sz2 == sc_sz) {
         Sc_Mat->CC_2_Num[temp_mat_dim_cc_2] = num2;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_cc_2 + temp_mat_dim_c_c + temp_mat_dim_cc_1], "D[DownUp(%d)]", tot_site/2);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_cc_2 + temp_mat_dim_c_c + temp_mat_dim_cc_1], "DownUp(r+1)");
         temp_mat_dim_cc_2++;
      }
   }
   
   
   Sc_Mat->dim_cc_1 = temp_mat_dim_cc_1;
   Sc_Mat->dim_cc_2 = temp_mat_dim_cc_2;
   Sc_Mat->dim_c_c  = temp_mat_dim_c_c;
   Sc_Mat->dim_tot  = temp_mat_dim_cc_1 + temp_mat_dim_cc_2 + temp_mat_dim_c_c;
   
   //Check Point
   if (temp_mat_dim_cc_1 != temp_mat_dim_cc_2) {
      printf("Error in MAKE_SC_MAT_BASIS_1DHUBBARD_VF\n");
      exit(1);
   }
   
}
