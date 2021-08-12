//Ascending order
#include "SML.h"

void QUICK_SORT_STABLE_INT3_SINT4(int *Base, int *Target, int *A1_Int, int *A2_Int, short *A1, short *A2, short *A3, short *A4, int left, int right) {
   
   if (right - left <= 1) {
      return;
   }
   
   int pivot_index = (left + right)/2;
   int pivot       = Target[pivot_index];
   int pivot_base  = Base[pivot_index];
   
   SWAP_INT(&Target[pivot_index], &Target[right - 1]);
   SWAP_INT(&Base[pivot_index]  , &Base[right - 1]  );
   SWAP_INT(&A1_Int[pivot_index], &A1_Int[right - 1] );
   SWAP_INT(&A2_Int[pivot_index], &A2_Int[right - 1] );
   SWAP_SINT(&A1[pivot_index]   , &A1[right - 1]    );
   SWAP_SINT(&A2[pivot_index]   , &A2[right - 1]    );
   SWAP_SINT(&A3[pivot_index]   , &A3[right - 1]    );
   SWAP_SINT(&A4[pivot_index]   , &A4[right - 1]    );
   
   int i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_INT(&Target[i], &Target[j]);
         SWAP_INT(&Base[i]  , &Base[j]  );
         SWAP_INT(&A1_Int[i], &A1_Int[j] );
         SWAP_INT(&A2_Int[i], &A2_Int[j] );
         SWAP_SINT(&A1[i]   , &A1[j]    );
         SWAP_SINT(&A2[i]   , &A2[j]    );
         SWAP_SINT(&A3[i]   , &A3[j]    );
         SWAP_SINT(&A4[i]   , &A4[j]    );
         i++;
      }
      else if (Target[j] == pivot) {
         if (Base[j] < pivot_base) {
            SWAP_INT(&Target[i], &Target[j]);
            SWAP_INT(&Base[i]  , &Base[j]  );
            SWAP_INT(&A1_Int[i], &A1_Int[j] );
            SWAP_INT(&A2_Int[i], &A2_Int[j] );
            SWAP_SINT(&A1[i]   , &A1[j]    );
            SWAP_SINT(&A2[i]   , &A2[j]    );
            SWAP_SINT(&A3[i]   , &A3[j]    );
            SWAP_SINT(&A4[i]   , &A4[j]    );
            i++;
         }
      }
   }
   
   SWAP_INT(&Target[i], &Target[right - 1]);
   SWAP_INT(&Base[i]  , &Base[right - 1]  );
   SWAP_INT(&A1_Int[i], &A1_Int[right - 1] );
   SWAP_INT(&A2_Int[i], &A2_Int[right - 1] );
   SWAP_SINT(&A1[i]   , &A1[right - 1]    );
   SWAP_SINT(&A2[i]   , &A2[right - 1]    );
   SWAP_SINT(&A3[i]   , &A3[right - 1]    );
   SWAP_SINT(&A4[i]   , &A4[right - 1]    );
   
   QUICK_SORT_STABLE_INT3_SINT4(Base, Target, A1_Int, A2_Int, A1, A2, A3, A4, left, i    );
   QUICK_SORT_STABLE_INT3_SINT4(Base, Target, A1_Int, A2_Int, A1, A2, A3, A4, i+1 , right);
   
   for (i = left; i < right; i++) {
      Base[i] = i;
   }
   
}
