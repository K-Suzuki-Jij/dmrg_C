#include <stdio.h>
#include <stdlib.h>
#include "SML.h"

CRS1 **GET_CRS2(int row, int dim, long max) {
   
   CRS1 **Matrix = malloc(sizeof(CRS1*)*row);
   
   int num;
   
   for (num = 0; num < row; num++) {
      Matrix[num] = GET_CRS1(dim, max);
   }
   
   return Matrix;
   
}
