//Ascending order
#include "SML.h"

void QUICK_SORT_LINT2_DOUBLE1(long *Target, long *A1, double *D2, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long pivot_index = (left + right)/2;
   long pivot       = Target[pivot_index];
   
   SWAP_LINT(&Target[pivot_index], &Target[right - 1]);
   SWAP_LINT(&A1[pivot_index]    , &A1[right - 1]    );
   SWAP_DOUBLE(&D2[pivot_index] , &D2[right - 1]    );
   
   long i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_LINT(&Target[i], &Target[j]);
         SWAP_LINT(&A1[i]    , &A1[j]    );
         SWAP_DOUBLE(&D2[i] , &D2[j]    );
         i++;
      }
   }
   
   SWAP_LINT(&Target[i], &Target[right - 1]);
   SWAP_LINT(&A1[i]    , &A1[right - 1]    );
   SWAP_DOUBLE(&D2[i] , &D2[right - 1]    );
   
   QUICK_SORT_LINT2_DOUBLE1(Target, A1, D2, left, i    );
   QUICK_SORT_LINT2_DOUBLE1(Target, A1, D2, i+1 , right);
   
}

