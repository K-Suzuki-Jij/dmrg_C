#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

double _Complex *GET_ARRAY_C_DOUBLE1(long num) {
   
   long i;
   
   double _Complex *Matrix = malloc(sizeof(double _Complex)*num);
   
   if (Matrix == NULL) {
      printf("Error in GET_ARRAY_C_DOUBLE1\n");
      printf("Need More Memory(num=%ld)\n", num);
      exit(1);
   }
   
   for (i = 0; i < num; i++) {
      Matrix[i] = 0;
   }
   
   return Matrix;
}
