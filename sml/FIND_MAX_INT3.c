#include <stdio.h>
#include <limits.h>

int FIND_MAX_INT3(int ***Array, long row1, long row2, long row3) {
   
   long i,j,k;
   long max = LONG_MIN;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         for (k = 0; k < row3; k++) {
            if (max < Array[i][j][k]) {
               max = Array[i][j][k];
            }
         }
      }
   }
   
   return (int)max;
   
}
