#include <stdio.h>
#include <limits.h>

long FIND_MAX_LINT1(long *Array, long dim) {
   
   long i;
   long max = LONG_MIN;
   
   for (i = 0; i < dim; i++) {
      if (max < Array[i]) {
         max = Array[i];
      }
   }
   
   return max;
   
}
