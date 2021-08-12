#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "SML.h"

void MATRIX_TRANSPOSE_CRS1(CRS1 *M, CRS1 *Out) {
   
   long i,j;
   long row_dim = M->row_dim;
   long col_dim = M->col_dim;
   long count1 = 0;
   long count2 = 0;
   int *Row_Count = GET_ARRAY_INT1(row_dim);
   FILE *file;
   
   //Bubble sort FIX ME:Replace marge sort
   for (j = 0; j < col_dim; j++) {
      count2 = 0;
      for (i = 0; i < row_dim; i++) {
         if ((M->Row[i] + Row_Count[i] < M->Row[i+1]) && (M->Col[M->Row[i] + Row_Count[i]] == j)) {
            Out->Val[count1] = M->Val[M->Row[i] + Row_Count[i]];
            Out->Col[count1] = (int)i;
            count1 = count1 + 1;
            count2 = count2 + 1;
            Row_Count[i] = Row_Count[i] + 1;
         }
      }
      
      if (j+1 >= Out->max_row) {
         printf("error in MATRIX_TRANSPOSE_CRS1\n");
         printf("j+1=%ld > out=%d\n", j+1, Out->max_row);
         mkdir("SML_out",0777);
         if ((file = fopen("./SML_out/error.txt","a+")) == NULL) {
            printf("Can't open file.\n");
            exit(1);
         }
         fprintf(file,"error in MATRIX_TRANSPOSE_CRS1\n");
         fclose(file);
         exit(1);
      }
      Out->Row[j+1] = Out->Row[j] + count2;
   }
   
   FREE_ARRAY_INT1(Row_Count);
   Out->row_dim = M->col_dim;
   Out->col_dim = M->row_dim;
   
}
