//Ascending order
#include "SML.h"

void QUICK_SORT_LINT1(long *Target, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long pivot_index = (left + right)/2;
   long pivot       = Target[pivot_index];
   
   SWAP_LINT(&Target[pivot_index], &Target[right - 1]);
   
   long i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_LINT(&Target[i], &Target[j]);
         i++;
      }
   }
   
   SWAP_LINT(&Target[i], &Target[right - 1]);
   
   QUICK_SORT_LINT1(Target, left, i    );
   QUICK_SORT_LINT1(Target, i+1 , right);
   
}

