#include <stdio.h>
#include <stdlib.h>
#include "SML.h"
#include "exact.h"

EXACT_A_BASIS **EXACT_GET_A_BASIS(int p_threads, int max_row) {
   
   EXACT_A_BASIS **A_Basis = malloc(sizeof(EXACT_A_BASIS*)*p_threads);
   int thread_num,i;
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      A_Basis[thread_num] = malloc(sizeof(EXACT_A_BASIS));
      A_Basis[thread_num]->Basis    = GET_ARRAY_LINT1(max_row);
      A_Basis[thread_num]->Check    = GET_ARRAY_LINT1(max_row);
      A_Basis[thread_num]->Val      = GET_ARRAY_DOUBLE1(max_row);
      A_Basis[thread_num]->max_row  = max_row;
      A_Basis[thread_num]->elem_num = 0;
      for (i = 0; i < max_row; i++) {
         A_Basis[thread_num]->Check[i] = -1;
      }
   }
   
   return A_Basis;
   
}


