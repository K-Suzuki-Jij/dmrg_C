#include <stdio.h>
#include <stdlib.h>
#include "SML.h"

CRS1 ***GET_CRS3(int row1, int row2, int dim, long max) {
   
   CRS1 ***Matrix = malloc(sizeof(CRS1**)*row1);
   
   int num;
   
   for (num = 0; num < row1; num++) {
      Matrix[num] = GET_CRS2(row2, dim, max);
   }
   
   return Matrix;
   
}
