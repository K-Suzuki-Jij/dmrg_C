#include <stdio.h>
#include <stdlib.h>

char ***GET_ARRAY_CHAR3(long row, long col, long col_2) {
   
   long i,j,k;
   
   char ***Matrix;
   Matrix = malloc(sizeof(char**)*row);
   if (Matrix == NULL) {
      printf("Error in GET_ARRAY_CHAR3\n");
      printf("Need More Memory(row=%ld)\n",row);
      exit(1);
   }
   for (i = 0; i < row; i++) {
      Matrix[i] = malloc(sizeof(char*)*col);
      if (Matrix[i] == NULL) {
         printf("Error in GET_ARRAY_CHAR3\n");
         printf("Need More Memory(col=%ld)\n", col);
         exit(1);
      }
      for (j = 0; j < col; j++) {
         Matrix[i][j] = malloc(sizeof(char)*col_2);
         if (Matrix[i][j] == NULL) {
            printf("Error in GET_ARRAY_CHAR3\n");
            printf("Need More Memory(col_2=%ld)\n", col_2);
            exit(1);
         }
      }
   }
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         for (k = 0; k < col_2; k++) {
            Matrix[i][j][k] = 0;
         }
      }
   }
   
   return Matrix;
   
}
