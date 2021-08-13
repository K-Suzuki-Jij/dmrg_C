#include <omp.h>
#include <math.h>
#include "dmrg.h"

void DMRG_V_M_LL_Q0(CRS1 *M_LL, double *V, double *Out_V, int dim, int p_threads, DMRG_BASIS *Dmrg_Basis) {
   
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
      for (iter = M_LL->Row[LL]; iter < M_LL->Row[LL + 1]; iter++) {
         col     = M_LL->Col[iter];
         inv_sup = (long)col*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_LL->Val[iter];
         }
      }
      Out_V[basis] = val;
   }

}
