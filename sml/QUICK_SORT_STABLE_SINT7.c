//Ascending order
#include "SML.h"

void QUICK_SORT_STABLE_SINT7(int *Base, short *Target, short *A1, short *A2, short *A3, short *A4, short *A5, short *A6, int left, int right) {
   
   if (right - left <= 1) {
      return;
   }
   
   int pivot_index = (left + right)/2;
   short pivot     = Target[pivot_index];
   int pivot_base  = Base[pivot_index];

   SWAP_INT(&Base[pivot_index]   , &Base[right - 1]  );
   SWAP_SINT(&Target[pivot_index], &Target[right - 1]);
   SWAP_SINT(&A1[pivot_index]    , &A1[right - 1]    );
   SWAP_SINT(&A2[pivot_index]    , &A2[right - 1]    );
   SWAP_SINT(&A3[pivot_index]    , &A3[right - 1]    );
   SWAP_SINT(&A4[pivot_index]    , &A4[right - 1]    );
   SWAP_SINT(&A5[pivot_index]    , &A5[right - 1]    );
   SWAP_SINT(&A6[pivot_index]    , &A6[right - 1]    );
   
   int i,j;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      if (Target[j] < pivot) {
         SWAP_INT(&Base[i]   , &Base[j]  );
         SWAP_SINT(&Target[i], &Target[j]);
         SWAP_SINT(&A1[i]    , &A1[j]    );
         SWAP_SINT(&A2[i]    , &A2[j]    );
         SWAP_SINT(&A3[i]    , &A3[j]    );
         SWAP_SINT(&A4[i]    , &A4[j]    );
         SWAP_SINT(&A5[i]    , &A5[j]    );
         SWAP_SINT(&A6[i]    , &A6[j]    );
         i++;
      }
      else if (Target[j] == pivot) {
         if (Base[j] < pivot_base) {
            SWAP_INT(&Base[i]   , &Base[j]  );
            SWAP_SINT(&Target[i], &Target[j]);
            SWAP_SINT(&A1[i]    , &A1[j]    );
            SWAP_SINT(&A2[i]    , &A2[j]    );
            SWAP_SINT(&A3[i]    , &A3[j]    );
            SWAP_SINT(&A4[i]    , &A4[j]    );
            SWAP_SINT(&A5[i]    , &A5[j]    );
            SWAP_SINT(&A6[i]    , &A6[j]    );
            i++;
         }
      }
   }

   SWAP_INT(&Base[i]   , &Base[right - 1]  );
   SWAP_SINT(&Target[i], &Target[right - 1]);
   SWAP_SINT(&A1[i]    , &A1[right - 1]    );
   SWAP_SINT(&A2[i]    , &A2[right - 1]    );
   SWAP_SINT(&A3[i]    , &A3[right - 1]    );
   SWAP_SINT(&A4[i]    , &A4[right - 1]    );
   SWAP_SINT(&A5[i]    , &A5[right - 1]    );
   SWAP_SINT(&A6[i]    , &A6[right - 1]    );
   
   QUICK_SORT_STABLE_SINT7(Base, Target, A1, A2, A3, A4, A5, A6, left, i    );
   QUICK_SORT_STABLE_SINT7(Base, Target, A1, A2, A3, A4, A5, A6, i+1 , right);
   
   for (i = left; i < right; i++) {
      Base[i] = i;
   }
   
}
