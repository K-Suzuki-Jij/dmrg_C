#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "SML.h"

double *GET_RAND_DOUBLE1(long num) {
   
   long i;
   
   double *Matrix = malloc(sizeof(double)*num);
   if (Matrix == NULL) {
      printf("Error in GET_RAND_DOUBLE1\n");
      printf("Need More Memory(num=%ld)\n", num);
      exit(1);
   }
   
   srand((unsigned int)time(NULL));
   
   for(i = 0;i < num;i++){
      Matrix[i] = rand()%10000000 - 5000000;
   }
   
   NORMALIZE(Matrix,num,1);
   
   return Matrix;
}
