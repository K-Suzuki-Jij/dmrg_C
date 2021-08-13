#include <stdlib.h>
#include "dmrg.h"

DMRG_A_BASIS **DMRG_GET_A_BASIS(int max_elem_num, int p_threads) {
   
   int thread_num;
   DMRG_A_BASIS **A_Basis = malloc(sizeof(DMRG_A_BASIS*)*p_threads);
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      A_Basis[thread_num] = malloc(sizeof(DMRG_A_BASIS));
      A_Basis[thread_num]->LL_LLLRRRRL  = GET_ARRAY_SINT1(max_elem_num);
      A_Basis[thread_num]->LR_LLLRRRRL  = GET_ARRAY_SINT1(max_elem_num);
      A_Basis[thread_num]->RR_LLLRRRRL  = GET_ARRAY_SINT1(max_elem_num);
      A_Basis[thread_num]->RL_LLLRRRRL  = GET_ARRAY_SINT1(max_elem_num);
      A_Basis[thread_num]->Val_LLLRRRRL = GET_ARRAY_DOUBLE1(max_elem_num);
      A_Basis[thread_num]->Inv_LLLRRRRL = GET_ARRAY_INT1(max_elem_num);
      A_Basis[thread_num]->elem_num = 0;
   }
   
   int basis;
   
#pragma omp parallel for private (basis) num_threads (p_threads)
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      for (basis = 0; basis < max_elem_num; basis++) {
         A_Basis[thread_num]->Inv_LLLRRRRL[basis] = -1;
      }
   }
   
   return A_Basis;
   
}
