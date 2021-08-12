#include <stdlib.h>
#include "SML.h"

CCS1 **GET_CCS2(int row, int dim, long max) {
   
   CCS1 **Matrix = malloc(sizeof(CCS1*)*row);
 
   int num;
   
   for (num = 0; num < row; num++) {
      Matrix[num] = GET_CCS1(dim, max);
   }
   
   return Matrix;
   
}
