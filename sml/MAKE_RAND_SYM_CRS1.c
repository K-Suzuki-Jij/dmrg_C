#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include "SML.h"

void MAKE_RAND_SYM_CRS1(CRS1 *M, long dim, long ele_num) {
   
   if (ele_num < dim || ele_num > dim*dim) {
      printf("Error in MAKE_RAND_SYM_CRS1\n");
      printf("ele_num=%ld,dim=%ld,dim*dim=%ld\n", ele_num, dim, dim*dim);
      exit(1);
   }
   
   if (M->max_val < ele_num || M->max_row < dim) {
      printf("Error in MAKE_RAND_SYM_CRS1\n");
      printf("M->max_val(%ld) < ele_num(%ld) or M->max_row(%d) < dim(%ld)\n", M->max_val, ele_num, M->max_row, dim);
      exit(1);
   }   
   
   long i,j,count,count2,flag;
   long diag_ele_num,i_temp;
   long *Off_Row_Ele_Num = GET_ARRAY_LINT1(dim);
   double rate;
   
   for (i = 0; i < M->max_val; i++) {
      M->Val[i] = 0;
      M->Col[i] = 0;
   }
   
   for (i = 0; i < M->max_row; i++) {
      M->Row[i] = 0;
   }
   
   M->row_dim = 0;
   M->col_dim = 0;

   
   if (ele_num%2 == 0) {
      do {
         diag_ele_num = rand()%(dim+1);
      } while (diag_ele_num%2 == 1);
   }
   else {
      do {
         diag_ele_num = rand()%(dim+1);
      } while (diag_ele_num%2 == 0);
   }

   i_temp = (ele_num - dim)/2;
   rate = (double)i_temp/(dim*(dim-1)/2);
   count = 0;
   for (i = 0; i < dim; i++) {
      Off_Row_Ele_Num[i] = (double)(dim - i - 1)*rate;
      count = count + Off_Row_Ele_Num[i];
   }
   
   flag = 0;
   
   for (i = 0; i < dim; i++) {
     
      if (2*count + diag_ele_num == ele_num || flag == 1) {
         break;
      }
      flag = 1;
      if (Off_Row_Ele_Num[i] < dim - i - 1) {
         Off_Row_Ele_Num[i] = Off_Row_Ele_Num[i] + 1;
         count = count + 1;
         flag = 0;
      }
   }

   if (2*count + diag_ele_num != ele_num) {
      diag_ele_num = ele_num - 2*count;
   }
   
   long off_ele_num = (ele_num - diag_ele_num)/2;
   int  *Temp_Col   = GET_ARRAY_INT1(off_ele_num);
   char *Col_Check  = GET_ARRAY_CHAR1(dim);

   count = 0;
   for (i = 0; i < dim; i++) {
      count2 = count;
      for (j = 0; j < Off_Row_Ele_Num[i]; j++) {
         Temp_Col[count] = (int)RANDOM(i+1,dim-1);
         if (Col_Check[Temp_Col[count]] == 1) {
            j = j - 1;
         }
         else {
            Col_Check[Temp_Col[count]] = 1;
            count = count + 1;
         }
      }
     
      for (j = count2; j < count; j++) {
         Col_Check[Temp_Col[j]] = 0;
      }
      if (count == off_ele_num) {
         break;
      }
   }

   
   int *Diag_Row = GET_ARRAY_INT1(diag_ele_num);
   for (i = 0; i < diag_ele_num; i++) {
      Diag_Row[i] = rand()%dim;
      if (Col_Check[Diag_Row[i]] == 1) {
         i = i - 1;
      }
      else {
         Col_Check[Diag_Row[i]] = 1;
      }
   }
   free(Col_Check);

   int *Row = GET_ARRAY_INT1(ele_num+1);

   count = 0;
   count2 = 0;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < Off_Row_Ele_Num[i]; j++) {
         M->Val[count] = rand()%10000000 - 5000000;
         M->Col[count] = Temp_Col[count2];
         Row[count]    = (int)i;
         M->Row[i+1]   = M->Row[i+1] + 1;
         count++;
         
         M->Val[count]        = M->Val[count-1];
         M->Col[count]        = (int)i;
         Row[count]           = Temp_Col[count2];
         M->Row[Row[count]+1] = M->Row[Row[count]+1] + 1;
         count++;
         count2++;
      }
   }
   
   free(Off_Row_Ele_Num);
   free(Temp_Col);

   for (i = 0; i < diag_ele_num; i++) {
      M->Val[count] = rand()%10000000 - 5000000;
      M->Col[count] = Diag_Row[i];
      Row[count] = Diag_Row[i];
      M->Row[Diag_Row[i]+1] = M->Row[Diag_Row[i]+1] + 1;
      count = count + 1;
   }

   free(Diag_Row);
   QUICK_SORT_INT2_DOUBLE1(Row, M->Col, M->Val, 0, ele_num);
   free(Row);
   
   for (i = 0; i < dim; i++) {
      M->Row[i+1] = M->Row[i+1] + M->Row[i];
   }
   
   M->row_dim = (int)dim;
   M->col_dim = (int)dim;
   
   NORMALIZE(M->Val, ele_num, 1);
   
   SORT_COLUMN_CRS1(M,1);

}
