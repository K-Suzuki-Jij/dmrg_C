#include "SML.h"
#include <stdlib.h>

void DMRG_QUICK_SORT_BASIS_Q1(short *Array_1, short *Array_2, short *Array_3, short *Array_4, short *Array_5, int left, int right) {
   
   if (right - left <= 1) {
      return;
   }
   
   int pivot_index = (left + right)/2;
   short pivot_1   = Array_1[pivot_index];
   short pivot_2   = Array_2[pivot_index];
   short pivot_3   = Array_3[pivot_index];
   short pivot_4   = Array_4[pivot_index];
   short pivot_5   = Array_5[pivot_index];
   
   SWAP_SINT(&Array_1[pivot_index], &Array_1[right - 1]);
   SWAP_SINT(&Array_2[pivot_index], &Array_2[right - 1]);
   SWAP_SINT(&Array_3[pivot_index], &Array_3[right - 1]);
   SWAP_SINT(&Array_4[pivot_index], &Array_4[right - 1]);
   SWAP_SINT(&Array_5[pivot_index], &Array_5[right - 1]);
   
   int i,j;
   int c1,c2,c3,c4,c5;
   int e2,e3,e4,e5;
   
   i = left;
   for (j = left; j < right - 1; j++) {
      
      e2 = (Array_2[j] == pivot_2);
      e3 = (Array_3[j] == pivot_3);
      e4 = (Array_4[j] == pivot_4);
      e5 = (Array_5[j] == pivot_5);

      c5 = (Array_5[j] < pivot_5);
      c4 = e5 && (Array_4[j] < pivot_4);
      c3 = e5 && e4 && (Array_3[j] < pivot_3);
      c2 = e5 && e4 && e3 && (Array_2[j] < pivot_2);
      c1 = e5 && e4 && e3 && e2 && (Array_1[j] < pivot_1);
      
      if (c1 || c2 || c3 || c4 || c5) {
         SWAP_SINT(&Array_1[i], &Array_1[j]);
         SWAP_SINT(&Array_2[i], &Array_2[j]);
         SWAP_SINT(&Array_3[i], &Array_3[j]);
         SWAP_SINT(&Array_4[i], &Array_4[j]);
         SWAP_SINT(&Array_5[i], &Array_5[j]);
         i++;
      }
   }
   
   SWAP_SINT(&Array_1[i], &Array_1[right - 1]);
   SWAP_SINT(&Array_2[i], &Array_2[right - 1]);
   SWAP_SINT(&Array_3[i], &Array_3[right - 1]);
   SWAP_SINT(&Array_4[i], &Array_4[right - 1]);
   SWAP_SINT(&Array_5[i], &Array_5[right - 1]);
   
   DMRG_QUICK_SORT_BASIS_Q1(Array_1, Array_2, Array_3, Array_4, Array_5, left, i    );
   DMRG_QUICK_SORT_BASIS_Q1(Array_1, Array_2, Array_3, Array_4, Array_5, i+1 , right);
   
}
