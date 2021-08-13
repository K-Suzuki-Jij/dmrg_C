#include <stdio.h>
#include <stdlib.h>
#include "SML.h"
#include "exact.h"

void EXACT_FREE_A_BASIS(EXACT_A_BASIS **A_Basis, int p_threads) {
   
   int thread_num;
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      FREE_ARRAY_LINT1(A_Basis[thread_num]->Basis);
      FREE_ARRAY_LINT1(A_Basis[thread_num]->Check);
      FREE_ARRAY_DOUBLE1(A_Basis[thread_num]->Val);
      free(A_Basis[thread_num]);
   }
   free(A_Basis);

}


