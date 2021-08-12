#include <stdio.h>
#include <stdlib.h>
#include "SML.h"

CRS1 ****GET_CRS4(int row1, int row2, int row3, int dim, long max) {
   
   CRS1 ****Matrix = malloc(sizeof(CRS1***)*row1);
   
   int num;
   
   for (num = 0; num < row1; num++) {
      Matrix[num] = GET_CRS3(row2, row3, dim, max);
   }
   
   return Matrix;
   
}
