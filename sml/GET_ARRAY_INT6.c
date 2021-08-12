#include <stdio.h>
#include <stdlib.h>

int ******GET_ARRAY_INT6(long row, long col, long col_2, long col_3, long col_4, long col_5) {
   
   
   long i,j,k,l,m,n;
   
   int ******Matrix;
   Matrix = malloc(sizeof(int*****)*row);
   if (Matrix == NULL) {
      printf("Error in GET_ARRAY_INT6\n");
      printf("Need More Memory(row=%ld)\n",row);
      exit(1);
   }
   for (i = 0; i < row; i++) {
      Matrix[i] = malloc(sizeof(int****)*col);
      if (Matrix[i] == NULL) {
         printf("Error in GET_ARRAY_INT6\n");
         printf("Need More Memory(col=%ld)\n",col);
         exit(1);
      }
      for (j = 0; j < col; j++) {
         Matrix[i][j] = malloc(sizeof(int***)*col_2);
         if (Matrix[i][j] == NULL) {
            printf("Error in GET_ARRAY_INT6\n");
            printf("Need More Memory(col_2=%ld)\n", col_2);
            exit(1);
         }
         for (k = 0; k < col_2; k++) {
            Matrix[i][j][k] = malloc(sizeof(int**)*col_3);
            if(Matrix[i][j][k] == NULL){
               printf("Error in GET_ARRAY_INT6\n");
               printf("Need More Memory(col_3=%ld)\n", col_3);
               exit(1);
            }
            for (l = 0; l < col_3; l++) {
               Matrix[i][j][k][l] = malloc(sizeof(int*)*col_4);
               if(Matrix[i][j][k][l] == NULL){
                  printf("Error in GET_ARRAY_INT6\n");
                  printf("Need More Memory(col_4=%ld)\n", col_4);
                  exit(1);
               }
               for (m = 0; m < col_4; m++) {
                  Matrix[i][j][k][l][m] = malloc(sizeof(int)*col_5);
                  if(Matrix[i][j][k][l][m] == NULL){
                     printf("Error in GET_ARRAY_INT6\n");
                     printf("Need More Memory(col_5=%ld)\n", col_5);
                     exit(1);
                  }
               }
            }
         }
      }
   }
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         for (k = 0; k < col_2; k++) {
            for (l = 0; l < col_3; l++) {
               for (m = 0; m < col_4; m++) {
                  for (n = 0; n < col_5; n++) {
                     Matrix[i][j][k][l][m][n] = 0;
                  }
               }
            }
         }
      }
   }
   
   return Matrix;
   
}
