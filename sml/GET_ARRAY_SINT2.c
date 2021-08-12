#include <stdio.h>
#include <stdlib.h>

short int **GET_ARRAY_SINT2(long row, long col) {
   
   long i,j;
   
   short int **Matrix;
   Matrix = malloc(sizeof(short int*)*row);
   if (Matrix == NULL) {
      printf("Error in GET_ARRAY_SINT2\n");
      printf("Need More Memory(row=%ld)\n", row);
      exit(1);
   }
   for (i = 0; i < row; i++) {
      Matrix[i] = malloc(sizeof(short int)*col);
      if (Matrix[i] == NULL) {
         printf("Error in GET_ARRAY_SINT2\n");
         printf("Need More Memory(col=%ld)\n", col);
         exit(1);
      }
   }
   
   for (i = 0; i < row; i++) {
      for (j = 0; j < col; j++) {
         Matrix[i][j] = 0;
      }
   }
   
   return Matrix;
   
}
