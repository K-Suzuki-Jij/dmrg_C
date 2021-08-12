#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include "SML.h"

void MAKE_RAND_CRS1(CRS1 *M, long dim, long ele_num) {

   if (dim <= 1) {
      printf("Error in MAKE_RAND_CRS1\n");
      printf("dim=%ld\n",dim);
      exit(1);
   }
   
   if (ele_num < dim || ele_num > dim*dim) {
      printf("Error in MAKE_RAND_CRS1\n");
      printf("ele_num=%ld,dim=%ld,dim*dim=%ld\n",ele_num, dim, dim*dim);
      exit(1);
   }
   
   if (M->max_val < ele_num || M->max_row < dim) {
      printf("Error in MAKE_SYM_CRS1\n");
      printf("M->max_val(%ld) < ele_num(%ld) or M->max_row(%d) < dim(%ld)\n", M->max_val, ele_num, M->max_row, dim);
      exit(1);
   }
      
   long i,j,count,count2;
   long *Row_Ele_Num = GET_ARRAY_LINT1(dim);
   
   for (i = 0; i < M->max_val; i++) {
      M->Val[i] = 0;
      M->Col[i] = 0;
   }
   
   for (i = 0; i < M->max_row; i++) {
      M->Row[i] = 0;
   }
   
   M->row_dim = 0;
   M->col_dim = 0;

   count = 0;
   for (i = 0; i < dim; i++) {
      Row_Ele_Num[i] = (double)ele_num/dim;
      count = count + Row_Ele_Num[i];
   }
   
   for (i = 0; i < dim; i++) {
      if (count == ele_num) {
         break;
      }
      Row_Ele_Num[i] = Row_Ele_Num[i] + 1;
      count = count + 1;
   }
   
   char *Col_Check = GET_ARRAY_CHAR1(dim);
   
   srand((unsigned int)time(NULL));
   count = 0;
   for (i = 0; i < dim; i++) {
      count2 = count;
      for (j = 0; j < Row_Ele_Num[i]; j++) {
         M->Col[count] = (int)rand()%(dim - 1);
         M->Val[count] = rand()%10000000 - 5000000;
         if (Col_Check[M->Col[count]] == 1) {
            j = j - 1;
         }
         else {
            Col_Check[M->Col[count]] = 1;
            count = count + 1;
         }
      }
      for (j = count2; j < count; j++) {
         Col_Check[M->Col[j]] = 0;
      }
      M->Row[i+1] = count;
      
      if (count == ele_num) {
         break;
      }
   }
   
   M->row_dim = (int)dim;
   M->col_dim = (int)dim;
   
   free(Col_Check);
   free(Row_Ele_Num);
   
   NORMALIZE(M->Val, ele_num, 1);
   SORT_COLUMN_CRS1(M, 1);
   
}
