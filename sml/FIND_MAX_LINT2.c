#include <stdio.h>
#include <limits.h>

long FIND_MAX_LINT2(long **Array, long row1, long row2) {
   
   long i,j;
   long max = LONG_MIN;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         if (max < Array[i][j]) {
            max = Array[i][j];
         }
      }
   }
   
   return max;
   
}
