#include <stdio.h>
#include <limits.h>

int FIND_MAX_INT4(int ****Array, long row1, long row2, long row3, long row4) {
   
   long i,j,k,l;
   long max = LONG_MIN;
   
   for (i = 0; i < row1; i++) {
      for (j = 0; j < row2; j++) {
         for (k = 0; k < row3; k++) {
            for (l = 0; l < row4; l++) {
               if (max < Array[i][j][k][l]) {
                  max = Array[i][j][k][l];
               }
            }
         }
      }
   }
   
   return (int)max;
   
}
