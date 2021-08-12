#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "SML.h"

void MATRIX_SUM_CRS1(CRS1 *M1, CRS1 *M2, CRS1 *Out) {
   
   if(M1->row_dim != M2->row_dim){
      printf("Error in MATRIX_SUM_CRS1\n");
      printf("M1_row dim(%d) does not equal to M2_row dim(%d)\n",M1->row_dim,M2->row_dim);
      exit(1);
   }
   if(M1->col_dim != M2->col_dim){
      printf("Error in MATRIX_SUM_CRS1\n");
      printf("M1_col dim(%d) does not equal to M2_col dim(%d)\n",M1->col_dim,M2->col_dim);
      exit(1);
   }

   FILE *file;
   long i,j;
   int row_dim = M1->row_dim;
   long m1_count,m2_count,total_count,count1,count2,check,row_count;
   double zero = pow(10,-15);
   
   total_count = 0;
   for (i = 0; i < row_dim; i++) {
      
      count1    = 0;
      count2    = 0;
      check     = 0;
      row_count = 0;
      m1_count  = M1->Row[i+1] - M1->Row[i];
      m2_count  = M2->Row[i+1] - M2->Row[i];
      
      if ((m1_count != 0)&&(m2_count == 0)) {
         for (j = M1->Row[i]; j < M1->Row[i+1]; j++) {
            if (total_count >= Out->max_val) {
               printf("Error in MATRIX_SUM_CRS1\n");
               printf("Need more Out->max_val=%ld\n", Out->max_val);
               mkdir("SML_out",0777);
               if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                  printf("Can't open file.\n");
                  exit(1);
               }
               fprintf(file, "Error in MATRIX_SUM_CRS1\n");
               fprintf(file, "Need more Out->max_val=%ld\n", Out->max_val);
               fclose(file);
               exit(1);
            }
            Out->Val[total_count] = M1->Val[j];
            Out->Col[total_count] = M1->Col[j];
            total_count = total_count + 1;
         }
         
         if (i+1 >= Out->max_row) {
            printf("Error in MATRIX_SUM_CRS1\n");
            printf("Need more Out->max_row=%ld\n", i+1);
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file, "Error in MATRIX_SUM_CRS1\n");
            fprintf(file, "Need more Out->max_row=%ld\n", i+1);
            fclose(file);
            exit(1);
         }
         Out->Row[i+1] = Out->Row[i] + M1->Row[i+1] - M1->Row[i];
      }
      
      else if ((m1_count == 0)&&(m2_count != 0)) {
         for (j = M2->Row[i]; j < M2->Row[i+1]; j++) {
            if (total_count >= Out->max_val) {
               printf("Error in MATRIX_SUM_CRS1\n");
               printf("Need more Out->max_val=%ld\n", Out->max_val);
               mkdir("SML_out",0777);
               if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                  printf("Can't open file.\n");
                  exit(1);
               }
               fprintf(file, "Error in MATRIX_SUM_CRS1\n");
               fprintf(file, "Need more Out->max_val=%ld\n", Out->max_val);
               fclose(file);
               exit(1);
            }
            Out->Val[total_count] = M2->Val[j];
            Out->Col[total_count] = M2->Col[j];
            total_count = total_count + 1;
         }
         
         if (i+1 >= Out->max_row) {
            printf("Error in MATRIX_SUM_CRS1\n");
            printf("Need more Out->max_row=%ld\n", i+1);
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file, "Error in MATRIX_SUM_CRS1\n");
            fprintf(file, "Need more Out->max_row=%ld\n", i+1);
            fclose(file);
            exit(1);
         }
         Out->Row[i+1] = Out->Row[i] + M2->Row[i+1] - M2->Row[i];
      }
      
      else if ((m1_count != 0)&&(m2_count != 0)) {
         for (j = 0; j < m1_count + m2_count; j++) {
            if (M1->Col[M1->Row[i] + count1] < M2->Col[M2->Row[i] + count2]) {
               if (total_count >= Out->max_val) {
                  printf("Error in MATRIX_SUM_CRS1\n");
                  printf("Need more Out->max_val=%ld\n", Out->max_val);
                  mkdir("SML_out",0777);
                  if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                     printf("Can't open file.\n");
                     exit(1);
                  }
                  fprintf(file, "Error in MATRIX_SUM_CRS1\n");
                  fprintf(file, "Need more Out->max_val=%ld\n", Out->max_val);
                  fclose(file);
                  exit(1);
               }
               Out->Val[total_count] = M1->Val[M1->Row[i] + count1];
               Out->Col[total_count] = M1->Col[M1->Row[i] + count1];
               total_count = total_count + 1;
               row_count = row_count + 1;
               count1 = count1 + 1;
               
               if (M1->Row[i] + count1 == M1->Row[i+1]) {
                  check = 1;
                  break;
               }
            }
            
            else if (M1->Col[M1->Row[i] + count1] == M2->Col[M2->Row[i] + count2]) {
               if (fabs(M1->Val[M1->Row[i] + count1] + M2->Val[M2->Row[i] + count2]) > zero) {
                  if (total_count >= Out->max_val) {
                     printf("Error in MATRIX_SUM_CRS1\n");
                     printf("Need more Out->max_val=%ld\n", Out->max_val);
                     mkdir("SML_out",0777);
                     if((file = fopen("./SML_out/error.txt","a+")) == NULL){
                        printf("Can't open file.\n");
                        exit(1);
                     }
                     fprintf(file,"Error in MATRIX_SUM_CRS1\n");
                     fprintf(file,"Need more Out->max_val=%ld\n", Out->max_val);
                     fclose(file);
                     exit(1);
                  }
                  Out->Val[total_count] = M1->Val[M1->Row[i] + count1] + M2->Val[M2->Row[i] + count2];
                  Out->Col[total_count] = M1->Col[M1->Row[i] + count1];
                  total_count = total_count + 1;
                  row_count = row_count + 1;
                  count1 = count1 + 1;
                  count2 = count2 + 1;
                  
                  if ((M1->Row[i] + count1 == M1->Row[i+1])&&(M2->Row[i] + count2 < M2->Row[i+1])) {
                     check = 1;
                     break;
                  }
                  else if ((M1->Row[i] + count1 < M1->Row[i+1])&&(M2->Row[i] + count2 == M2->Row[i+1])) {
                     check = 2;
                     break;
                  }
                  else if ((M1->Row[i] + count1 == M1->Row[i+1])&&(M2->Row[i] + count2 == M2->Row[i+1])) {
                     check = 0;
                     break;
                  }
               }
               
               else {
                  count1 = count1 + 1;
                  count2 = count2 + 1;
                  
                  if ((M1->Row[i] + count1 == M1->Row[i+1])&&(M2->Row[i] + count2 < M2->Row[i+1])) {
                     check = 1;
                     break;
                  }
                  else if ((M1->Row[i] + count1 < M1->Row[i+1])&&(M2->Row[i] + count2 == M2->Row[i+1])) {
                     check = 2;
                     break;
                  }
                  else if ((M1->Row[i] + count1 == M1->Row[i+1])&&(M2->Row[i] + count2 == M2->Row[i+1])) {
                     check = 0;
                     break;
                  }
               }
            }
            
            else {
               if (total_count >= Out->max_val) {
                  printf("Error in MATRIX_SUM_CRS1\n");
                  printf("Need more Out->max_val=%ld\n", Out->max_val);
                  mkdir("SML_out",0777);
                  if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                     printf("Can't open file.\n");
                     exit(1);
                  }
                  fprintf(file,"Error in MATRIX_SUM_CRS1\n");
                  fprintf(file,"Need more Out->max_val=%ld\n", Out->max_val);
                  fclose(file);
                  exit(1);
               }
               Out->Val[total_count] = M2->Val[M2->Row[i] + count2];
               Out->Col[total_count] = M2->Col[M2->Row[i] + count2];
               total_count = total_count + 1;
               row_count = row_count + 1;
               count2 = count2 + 1;
               
               if (M2->Row[i] + count2 == M2->Row[i+1]) {
                  check = 2;
                  break;
               }
            }
         }
         
         if (check == 1) {
            for (j = M2->Row[i] + count2; j < M2->Row[i+1]; j++) {
               if (total_count >= Out->max_val) {
                  printf("Error in MATRIX_SUM_CRS1\n");
                  printf("Need more Out->max_val=%ld\n",Out->max_val);
                  mkdir("SML_out",0777);
                  if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                     printf("Can't open file.\n");
                     exit(1);
                  }
                  fprintf(file, "Error in MATRIX_SUM_CRS1\n");
                  fprintf(file, "Need more Out->max_val=%ld\n", Out->max_val);
                  fclose(file);
                  exit(1);
               }
               Out->Val[total_count] = M2->Val[j];
               Out->Col[total_count] = M2->Col[j];
               total_count = total_count + 1;
               row_count = row_count + 1;
            }
         }
         
         else if (check == 2) {
            for (j = M1->Row[i] + count1; j < M1->Row[i+1]; j++) {
               if (total_count >= Out->max_val) {
                  printf("Error in MATRIX_SUM_CRS1\n");
                  printf("Need more Out->max_val=%ld\n", Out->max_val);
                  mkdir("SML_out",0777);
                  if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
                     printf("Can't open file.\n");
                     exit(1);
                  }
                  fprintf(file,"Error in MATRIX_SUM_CRS1\n");
                  fprintf(file,"Need more Out->max_val=%ld\n", Out->max_val);
                  fclose(file);
                  exit(1);
               }
               Out->Val[total_count] = M1->Val[j];
               Out->Col[total_count] = M1->Col[j];
               total_count = total_count + 1;
               row_count = row_count + 1;
            }
         }
         
         if (i+1 >= Out->max_row) {
            printf("Error in MATRIX_SUM_CRS1\n");
            printf("Need more Out->max_row=%ld\n", i+1);
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file,"Error in MATRIX_SUM_CRS1\n");
            fprintf(file,"Need more Out->max_row=%ld\n", i+1);
            fclose(file);
            exit(1);
         }
         Out->Row[i+1] = Out->Row[i] + row_count;
      }
      
      else {
         if (i+1 >= Out->max_row) {
            printf("Error in MATRIX_SUM_CRS1\n");
            printf("Need more Out->max_row=%ld\n", i+1);
            mkdir("SML_out",0777);
            if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
               printf("Can't open file.\n");
               exit(1);
            }
            fprintf(file,"Error in MATRIX_SUM_CRS1\n");
            fprintf(file,"Need more Out->max_row=%ld\n", i+1);
            fclose(file);
            exit(1);
         }
         Out->Row[i+1] = Out->Row[i];
      }
   }
   
   if (row_dim >= Out->max_row) {
      printf("Error in MATRIX_SUM_CRS1\n");
      printf("Need more Out->max_row=%ld\n", i+1);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"Error in MATRIX_SUM_CRS1\n");
      fprintf(file,"Need more Out->max_row=%ld\n", i+1);
      fclose(file);
      exit(1);
   }
   
   Out->Row[row_dim] = total_count;
   Out->row_dim = row_dim;
   Out->col_dim = M1->col_dim;
   
   if (total_count > Out->max_val) {
      printf("Error in MATRIX_SUM_CRS1\n");
      printf("Need more Out->max_val=%ld\n", Out->max_val);
      mkdir("SML_out",0777);
      if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
         printf("Can't open file.\n");
         exit(1);
      }
      fprintf(file,"Error in MATRIX_SUM_CRS1\n");
      fprintf(file,"Need more Out->max_val=%ld\n", Out->max_val);
      fclose(file);
      exit(1);
   }
   
}
