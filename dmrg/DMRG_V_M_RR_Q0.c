#include <omp.h>
#include <math.h>
#include "dmrg.h"

void DMRG_V_M_RR_Q0(CRS1 *M_RR, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis) {
   
   int dim_RR     = Dmrg_Basis->dim_RR;
   int dim_onsite = Dmrg_Basis->dim_onsite;

   int basis,LL,LR,RR,RL,col,inv;
   long iter,inv_sup;
   double val;
   
#pragma omp parallel for private (LL,LR,RR,RL,val,iter,col,inv,inv_sup) num_threads (p_threads)
   for (basis = 0; basis < dim; basis++) {
      LL = Dmrg_Basis->LL_LLLRRRRL[basis];
      LR = Dmrg_Basis->LR_LLLRRRRL[basis];
      RR = Dmrg_Basis->RR_LLLRRRRL[basis];
      RL = Dmrg_Basis->RL_LLLRRRRL[basis];
      val = 0;
      for (iter = M_RR->Row[RR]; iter < M_RR->Row[RR + 1]; iter++) {
         col     = M_RR->Col[iter];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + col*dim_onsite + RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_RR->Val[iter];
         }
      }
      Out_V[basis] = val;
   }
   
}
