#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

void MAKE_SC_MAT_BASIS_1DKLM_TVF(int sc_p, SC_MAT_1DKLM_TVF *Sc_Mat, int tot_site) {
   
   ///////////////////////////////////////
   // # <->  [Cherge   ] -- (N,  P)
   // 0 <->  [         ] -- (0,  0)
   // 1 <->  [even     ] -- (1,  0)
   // 2 <->  [odd      ] -- (1,  1)
   // 3 <->  [even, odd] -- (2,  1)
   ///////////////////////////////////////
   
   if (sc_p < 0 || sc_p > 1) {
      printf("MAKE_SC_MAT_BASIS_1DKLM_TVF\n");
      printf("sc_p=%d\n",sc_p);
      exit(1);
   }

   int dim_cc_onsite  = 1;
   int dim_c_onsite   = 2;
   int num1,num2,comp_p1,comp_p2,temp_mat_dim_cc_1,temp_mat_dim_cc_2,temp_mat_dim_c_c;
   char C_Name1[100], C_Name2[100];
   
   //Onsite SC correlations Even*Odd(i)
   temp_mat_dim_cc_1 = 0;
   for (num1 = 0; num1 < dim_cc_onsite; num1++) {
      comp_p1 = 1;
      if (comp_p1 == sc_p) {
         Sc_Mat->CC_1_Num[temp_mat_dim_cc_1] = num1;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_cc_1], "D[OddEven(%d)]", tot_site/2);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_cc_1], "OddEven(r)");
         temp_mat_dim_cc_1++;
      }
      Sc_Mat->CC_Parity[num1] = comp_p1;
   }
   
   //Intersite SC correlations C_1(i)C_2(i+1)
   temp_mat_dim_c_c = 0;
   for (num1 = 0; num1 < dim_c_onsite; num1++) {
    
      if (num1 == 0) {
         comp_p1 = 0;
         sprintf(C_Name1, "Even");
      }
      else if (num1 == 1) {
         comp_p1 = 1;
         sprintf(C_Name1, "Odd");
      }
      else {
         printf("Error in MAKE_SC_MAT_BASIS_1DKLM_TVF\n");
         printf("num1=%d\n",num1);
         exit(1);
      }
      
      for (num2 = 0; num2 < dim_c_onsite; num2++) {
         if (num2 == 0) {
            comp_p2 = 0;
            sprintf(C_Name2, "Even");
         }
         else if (num2 == 1) {
            comp_p2 = 1;
            sprintf(C_Name2, "Odd");
         }
         else {
            printf("Error in MAKE_SC_MAT_BASIS_1DKLM_VF\n");
            printf("num2=%d\n",num2);
            exit(1);
         }
         if ((comp_p1 + comp_p2)%2 == sc_p) {
            Sc_Mat->C_Num1[temp_mat_dim_c_c] = num1;
            Sc_Mat->C_Num2[temp_mat_dim_c_c] = num2;
            sprintf(Sc_Mat->Row_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1], "D[%s(%d)%s(%d)]", C_Name1, tot_site/2, C_Name2, tot_site/2 - 1);
            sprintf(Sc_Mat->Col_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1], "%s(r)%s(r+1)", C_Name1, C_Name2);
            temp_mat_dim_c_c++;
         }
      }
      Sc_Mat->C_Parity[num1] = comp_p1;
   }
   
   //Onsite SC correlations Even*Odd(i+1)
   temp_mat_dim_cc_2 = 0;
   for (num2 = 0; num2 < dim_cc_onsite; num2++) {
     comp_p2 = 1;
      if (comp_p2 == sc_p) {
         Sc_Mat->CC_2_Num[temp_mat_dim_cc_2] = num2;
         sprintf(Sc_Mat->Row_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1 + temp_mat_dim_cc_2], "D[OddEven(%d)]", tot_site/2 - 1);
         sprintf(Sc_Mat->Col_Name[temp_mat_dim_c_c + temp_mat_dim_cc_1 + temp_mat_dim_cc_2], "OddEven(r+1)");
         temp_mat_dim_cc_2++;
      }
   }
   
   Sc_Mat->dim_cc_1 = temp_mat_dim_cc_1;
   Sc_Mat->dim_cc_2 = temp_mat_dim_cc_2;
   Sc_Mat->dim_c_c  = temp_mat_dim_c_c;
   Sc_Mat->dim_tot  = temp_mat_dim_cc_1 + temp_mat_dim_cc_2 + temp_mat_dim_c_c;
   
   //Check Point
   if (temp_mat_dim_cc_1 != temp_mat_dim_cc_2) {
      printf("Error in MAKE_SC_MAT_BASIS_1DKLM_TVF\n");
      printf("temp_mat_dim_cc_1=%d\n",temp_mat_dim_cc_1);
      printf("temp_mat_dim_cc_2=%d\n",temp_mat_dim_cc_2);
      exit(1);
   }
   
}
