//Ascending order
#include "SML.h"

void QUICK_SORT_INT2_DOUBLE1(int *Target, int *A1, double *D2, long left, long right) {
   
   if (right - left <= 1) {
      return;
   }
   
   long pivot_index = (left + right)/2;
   int pivot       = Target[pivot_index];
   
   SWAP_INT(&Target[pivot_index], &Target[right - 1]);
   SWAP_INT(&A1[pivot_index]    , &A1[right - 1]    );
   SWAP_DOUBLE(&D2[pivot_index] , &D2[right - 1]    );

   long i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_INT(&Target[i], &Target[j]);
         SWAP_INT(&A1[i]    , &A1[j]    );
         SWAP_DOUBLE(&D2[i] , &D2[j]    );
         i++;
      }
   }
   
   SWAP_INT(&Target[i], &Target[right - 1]);
   SWAP_INT(&A1[i]    , &A1[right - 1]    );
   SWAP_DOUBLE(&D2[i] , &D2[right - 1]    );

   QUICK_SORT_INT2_DOUBLE1(Target, A1, D2, left, i    );
   QUICK_SORT_INT2_DOUBLE1(Target, A1, D2, i+1 , right);

}

