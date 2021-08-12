#include <stdio.h>
#include <limits.h>

int FIND_MAX_INT1(int *Array, long dim) {
   
   long i;
   int max = INT_MIN;
   
   for (i = 0; i < dim; i++) {
      if (max < Array[i]) {
         max = Array[i];
      }
   }
   
   return max;
   
}
