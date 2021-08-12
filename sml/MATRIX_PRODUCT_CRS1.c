#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "SML.h"

void MATRIX_PRODUCT_CRS1(CRS1 *M1, CRS1 *M2, CRS1 *Out) {
  
   if (M1->col_dim != M2->row_dim) {
      printf("error in MATRIX_PRODUCT_CRS1\n");
      printf("Matrix product can't be defined\n");
      printf("M1_col=%d,M2_row=%d\n",M1->col_dim,M2->row_dim);
      exit(1);
   }
   
   FILE *file;
   long i,j,k;
   int m1_row_dim = M1->row_dim;
   int m1_col_dim = M1->col_dim;
   int m2_col_dim = M2->col_dim;
   int total_count,row_count;
   double zero = pow(10,-15);
   double *Temp_1 = GET_ARRAY_DOUBLE1(m1_col_dim);
   double *Temp_2 = GET_ARRAY_DOUBLE1(m2_col_dim);
   
   total_count = 0;
   for (i = 0; i < m1_row_dim; i++) {
      
      for (j = M1->Row[i]; j < M1->Row[i+1]; j++) {
         Temp_1[M1->Col[j]] = M1->Val[j];
      }
      
      for (j = 0; j < m1_col_dim; j++) {
         if (fabs(Temp_1[j]) > zero) {
            for (k = M2->Row[j]; k < M2->Row[j+1]; k++) {
               Temp_2[M2->Col[k]] = Temp_2[M2->Col[k]] + Temp_1[j]*M2->Val[k];
            }
         }
      }
      
      row_count = 0;
      for (j = 0; j < m2_col_dim; j++) {
         if (fabs(Temp_2[j]) > zero) {
            Out->Val[total_count] = Temp_2[j];
            Out->Col[total_count] = (int)j;
            total_count = total_count + 1;
            row_count   = row_count + 1;
            if (total_count > Out->max_val) {
               printf("Error in MATRIX_PRODUCT_CRS1\n");
               printf("need more Out->max_val=%ld\n", Out->max_val);
               mkdir("SML_out",0777);
               if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                  printf("Can't open file.\n");
                  exit(1);
               }
               fprintf(file,"Error in MATRIX_PRODUCT_CRS1\n");
               fprintf(file,"need more Out->max_val=%ld\n", Out->max_val);
               fclose(file);
               exit(1);
            }
         }
      }
      
      for (j = 0; j < m2_col_dim; j++) {
         Temp_2[j] = 0;
      }
      
      for (j = M1->Row[i]; j < M1->Row[i+1]; j++) {
         Temp_1[M1->Col[j]] = 0;
      }
      
      if (i+1 >= Out->max_row) {
         printf("error in MATRIX_PRODUCT_CRS1\n");
         printf("need more Out->max_row=%d\n", Out->max_row);
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file,"Error in MATRIX_PRODUCT_CRS1\n");
         fprintf(file,"need more Out->max_row=%d\n", Out->max_row);
         fclose(file);
         exit(1);
      }
      Out->Row[i+1] = Out->Row[i] + row_count;
   }
   
   FREE_ARRAY_DOUBLE1(Temp_1);
   FREE_ARRAY_DOUBLE1(Temp_2);

   Out->row_dim = m1_row_dim;
   Out->col_dim = m2_col_dim;   
}

















