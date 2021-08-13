#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "dmrg.h"

void DMRG_MATRIX_VECTOR_PRODUCT_OBC(CRS1 *M_LLLR, CRS1 *M_LRRL, CRS1 *M_LRRL_Sign, CRS1 *M_RRRL, double *V, double *Out_V, int *Ele_RR, int p_threads, DMRG_BASIS *Dmrg_Basis) {
   
   int dim_RR     = Dmrg_Basis->dim_RR;
   int dim_onsite = Dmrg_Basis->dim_onsite;
   int basis,LLLR,LRRL,RRRL,col,inv,sign;
   int LL,LR,RR,RL;
   int A_LL,A_LR,A_RR,A_RL;
   long iter,inv_sup;
   double val;
   int dim = Dmrg_Basis->dim_LLLRRRRL;
   
#pragma omp parallel for private (LL,LR,RR,RL,LLLR,LRRL,RRRL,iter,col,A_LL,A_LR,A_RR,A_RL,sign,inv,val,inv_sup) schedule(guided) num_threads (p_threads)
   for (basis = 0; basis < dim; basis++) {
      LL   = Dmrg_Basis->LL_LLLRRRRL[basis];
      LR   = Dmrg_Basis->LR_LLLRRRRL[basis];
      RR   = Dmrg_Basis->RR_LLLRRRRL[basis];
      RL   = Dmrg_Basis->RL_LLLRRRRL[basis];
      LLLR = Dmrg_Basis->Inv_LLLR[LL][LR];
      LRRL = Dmrg_Basis->Inv_LRRL[LR][RL];
      RRRL = Dmrg_Basis->Inv_RRRL[RR][RL];
      val  = 0;

      //Ham_LLLR
      for (iter = M_LLLR->Row[LLLR]; iter < M_LLLR->Row[LLLR + 1]; iter++) {
         col     = M_LLLR->Col[iter];
         A_LL    = Dmrg_Basis->LL_LLLR[col];
         A_LR    = Dmrg_Basis->LR_LLLR[col];
         inv_sup = (long)A_LL*dim_onsite*dim_RR*dim_onsite + (long)A_LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_LLLR->Val[iter];
         }
      }
      
      //Ham_LRRL
      for (iter = M_LRRL->Row[LRRL]; iter < M_LRRL->Row[LRRL + 1]; iter++) {
         col     = M_LRRL->Col[iter];
         A_LR    = Dmrg_Basis->LR_LRRL[col];
         A_RL    = Dmrg_Basis->RL_LRRL[col];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)A_LR*dim_RR*dim_onsite + RR*dim_onsite + A_RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_LRRL->Val[iter];
         }
      }
      
      if (Ele_RR[RR]%2 == 1) {
         sign = -1;
      }
      else {
         sign = 1;
      }
      
      //Ham_LRRL_Sign
      for (iter = M_LRRL_Sign->Row[LRRL]; iter < M_LRRL_Sign->Row[LRRL + 1]; iter++) {
         col     = M_LRRL_Sign->Col[iter];
         A_LR    = Dmrg_Basis->LR_LRRL[col];
         A_RL    = Dmrg_Basis->RL_LRRL[col];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)A_LR*dim_RR*dim_onsite + RR*dim_onsite + A_RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_LRRL_Sign->Val[iter]*sign;
         }
      }
      
      //Ham_RRRL
      for (iter = M_RRRL->Row[RRRL]; iter < M_RRRL->Row[RRRL + 1]; iter++) {
         col     = M_RRRL->Col[iter];
         A_RR    = Dmrg_Basis->RR_RRRL[col];
         A_RL    = Dmrg_Basis->RL_RRRL[col];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + A_RR*dim_onsite + A_RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + V[inv]*M_RRRL->Val[iter];
         }
      }
      Out_V[basis] = val;

   }
   
}
