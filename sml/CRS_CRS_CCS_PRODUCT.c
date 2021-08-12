#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"

void CRS_CRS_CCS_PRODUCT(CRS1 *M1, CRS1 *M2, CCS1 *M3, CRS1 *Out, CCS1 *Work, double *Work_Vec) {
   
   //Check Point
   if (M1->col_dim != M2->row_dim || M2->col_dim != M3->row_dim) {
      printf("error in CRS_CRS_CCS_PRODUCT\n");
      printf("M1_col = %d, M2_row = %d, M2_col = %d, M3_col = %d\n", M1->col_dim, M2->row_dim, M2->col_dim, M3->row_dim);
      exit(1);
   }
   
   //Check Point
   if (Work->max_col <= M3->col_dim) {
      printf("error in CRS_CRS_CCS_PRODUCT\n");
      printf("Need more Work->max_col = %d <= %d\n", Work->max_col, M3->col_dim);
      exit(1);
   }
   
   //Check Point
   if (Out->max_row <= M1->row_dim) {
      printf("error in CRS_CRS_CCS_PRODUCT\n");
      printf("Need more Out->max_row = %d <= %d\n", Out->max_row, M1->row_dim);
      exit(1);
   }
   
   long elem_num = 0;
   int row_m1 = M1->row_dim;
   int row_m2 = M2->row_dim;
   int row_m3 = M3->row_dim;
   int col_m3 = M3->col_dim;
   long iter1,iter2,iter3;
   double zero = pow(10,-15);
   double val;
   
   for (iter1 = 0; iter1 < row_m3; iter1++) {
      Work_Vec[iter1] = 0;
   }
   
   //M2*M3 --> Out in CCS
   for (iter1 = 0; iter1 < col_m3; iter1++) {
      
      for (iter2 = M3->Col[iter1]; iter2 < M3->Col[iter1 + 1]; iter2++) {
         Work_Vec[M3->Row[iter2]] = M3->Val[iter2];
      }
      
      for (iter2 = 0; iter2 < row_m2; iter2++) {
         val = 0;
         for (iter3 = M2->Row[iter2]; iter3 < M2->Row[iter2 + 1]; iter3++) {
            val = val + Work_Vec[M2->Col[iter3]]*M2->Val[iter3];
         }
         
         //Check Point
         if (elem_num >= Work->max_val) {
            printf("error in CRS_CRS_CCS_PRODUCT\n");
            printf("Need more Work->max_val = %ld\n",elem_num);
            exit(1);
         }
         
         if (fabs(val) > zero) {
            Work->Val[elem_num] = val;
            Work->Row[elem_num] = (int)iter2;
            elem_num++;
         }
      }
      
      for (iter2 = M3->Col[iter1]; iter2 < M3->Col[iter1 + 1]; iter2++) {
         Work_Vec[M3->Row[iter2]] = 0.0;
      }
      
      Work->Col[iter1 + 1] = elem_num;
      
   }
   
   Work->row_dim = row_m2;
   Work->col_dim = col_m3;
   
   
   for (iter1 = 0; iter1 < row_m2; iter1++) {
      Work_Vec[iter1] = 0;
   }
   
   elem_num = 0;
   
   //M1*(M2*M3) -->Out in CRS
   for (iter1 = 0; iter1 < row_m1; iter1++) {
      
      for (iter2 = M1->Row[iter1]; iter2 < M1->Row[iter1 + 1]; iter2++) {
         Work_Vec[M1->Col[iter2]] = M1->Val[iter2];
      }
      
      for (iter2 = 0; iter2 < col_m3; iter2++) {
         val = 0;
         for (iter3 = Work->Col[iter2]; iter3 < Work->Col[iter2 + 1]; iter3++) {
            val = val + Work_Vec[Work->Row[iter3]]*Work->Val[iter3];
         }
         
         //Check Point
         if (elem_num >= Out->max_val) {
            printf("error in CRS_CRS_CCS_PRODUCT\n");
            printf("Need more Out->max_val = %ld\n",elem_num);
            exit(1);
         }
         
         if (fabs(val) > zero) {
            Out->Val[elem_num] = val;
            Out->Col[elem_num] = (int)iter2;
            elem_num++;
         }
      }
      
      for (iter2 = M1->Row[iter1]; iter2 < M1->Row[iter1 + 1]; iter2++) {
         Work_Vec[M1->Col[iter2]] = 0.0;
      }
      
      Out->Row[iter1 + 1] = elem_num;
   }
   
   Out->row_dim = row_m1;
   Out->col_dim = col_m3;
   
}
