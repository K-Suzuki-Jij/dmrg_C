#include "dmrg.h"
#include <stdlib.h>

void DMRG_RE_ALLOCATE_INV_LLLRRRRL(DMRG_BASIS *Basis, int dim_LL, int dim_RR, int dim_onsite, int p_threads) {
   
   int dim_LLLRRRRL = Basis->dim_LLLRRRRL;
   long dim_whole   = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   int basis,LL,LR,RR,RL;
   
   Basis->Inv_LLLRRRRL = GET_ARRAY_INT1(dim_whole);
   
#pragma omp parallel for num_threads (p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }
   
   
#pragma omp parallel for private (LL,LR,RL,RR,inv_sup) num_threads (p_threads)
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      LL = Basis->LL_LLLRRRRL[basis];
      LR = Basis->LR_LLLRRRRL[basis];
      RL = Basis->RL_LLLRRRRL[basis];
      RR = Basis->RR_LLLRRRRL[basis];
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }

}
