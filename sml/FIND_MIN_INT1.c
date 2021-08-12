#include <stdio.h>
#include <limits.h>

int FIND_MIN_INT1(int *Array, int dim) {
   
   long i;
   int min = INT_MAX;
   
   for (i = 0; i < dim; i++) {
      if (min > Array[i]) {
         min = Array[i];
      }
   }
   
   return min;
   
}
