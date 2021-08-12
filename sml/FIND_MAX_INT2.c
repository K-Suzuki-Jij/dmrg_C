#include <stdio.h>
#include <limits.h>

int FIND_MAX_INT2(int **Array, long row1, long row2) {
   
   long i,j;
   int max = INT_MIN;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         if (max < Array[i][j]) {
            max = Array[i][j];
         }
      }
   }
   
   return max;
   
}
