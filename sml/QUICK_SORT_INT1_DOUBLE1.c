//Ascending order
#include "SML.h"

void QUICK_SORT_INT1_DOUBLE1(int *Target, double *D1, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long pivot_index = (left + right)/2;
   int pivot       = Target[pivot_index];
   
   SWAP_INT(&Target[pivot_index], &Target[right - 1]);
   SWAP_DOUBLE(&D1[pivot_index] , &D1[right - 1]    );
  
   long i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_INT(&Target[i], &Target[j]);
         SWAP_DOUBLE(&D1[i] , &D1[j]    );
         i++;
      }
   }
   
   SWAP_INT(&Target[i], &Target[right - 1]);
   SWAP_DOUBLE(&D1[i] , &D1[right - 1]    );

   QUICK_SORT_INT1_DOUBLE1(Target, D1, left, i    );
   QUICK_SORT_INT1_DOUBLE1(Target, D1, i+1 , right);

}

