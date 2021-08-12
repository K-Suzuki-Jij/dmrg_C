#include <stdio.h>
#include <stdlib.h>

double *******GET_ARRAY_DOUBLE7(long row, long col, long col_2, long col_3, long col_4, long col_5, long col_6) {
   
   long i,j,k,l,m,n,o;
   
   double *******Matrix;
   Matrix = malloc(sizeof(double******)*row);
   if (Matrix == NULL) {
      printf("Error in GET_ARRAY_DOUBLE7\n");
      printf("Need More Memory(row=%ld)\n", row);
      exit(1);
   }
   for (i = 0; i < row; i++) {
      Matrix[i] = malloc(sizeof(double*****)*col);
      if (Matrix[i] == NULL) {
         printf("Error in GET_ARRAY_DOUBLE7\n");
         printf("Need More Memory(col=%ld)\n", col);
         exit(1);
      }
      for (j = 0; j < col; j++) {
         Matrix[i][j] = malloc(sizeof(double****)*col_2);
         if (Matrix[i][j] == NULL) {
            printf("Error in GET_ARRAY_DOUBLE7\n");
            printf("Need More Memory(col_2=%ld)\n", col_2);
            exit(1);
         }
         for (k = 0; k < col_2; k++) {
            Matrix[i][j][k] = malloc(sizeof(double***)*col_3);
            if(Matrix[i][j][k] == NULL){
               printf("Error in GET_ARRAY_DOUBLE7\n");
               printf("Need More Memory(col_3=%ld)\n", col_3);;
               exit(1);
            }
            for (l = 0; l < col_3; l++) {
               Matrix[i][j][k][l] = malloc(sizeof(double**)*col_4);
               if(Matrix[i][j][k][l] == NULL){
                  printf("Error in GET_ARRAY_DOUBLE7\n");
                  printf("Need More Memory(col_4=%ld)\n", col_4);
                  exit(1);
               }
               for (m = 0; m < col_4; m++) {
                  Matrix[i][j][k][l][m] = malloc(sizeof(double*)*col_5);
                  if(Matrix[i][j][k][l][m] == NULL){
                     printf("Error in GET_ARRAY_DOUBLE7\n");
                     printf("Need More Memory(col_5=%ld)\n", col_5);
                     exit(1);
                  }
                  for (n = 0; n < col_5; n++) {
                     Matrix[i][j][k][l][m][n] = malloc(sizeof(double)*col_6);
                     if(Matrix[i][j][k][l][m][n] == NULL){
                        printf("Error in GET_ARRAY_DOUBLE7\n");
                        printf("Need More Memory(col_6=%ld)\n", col_6);
                        exit(1);
                     }
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
                     for (o = 0; o < col_6; o++) {
                        Matrix[i][j][k][l][m][n][o] = 0;
                     }
                  }
               }
            }
         }
      }
   }
   
   return Matrix;
   
}
