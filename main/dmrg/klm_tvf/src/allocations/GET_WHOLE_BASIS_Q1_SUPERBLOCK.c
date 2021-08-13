//
//  GET_WHOLE_BASIS_Q1_SUPERBLOCK.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/17.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_WHOLE_BASIS_Q1 *GET_WHOLE_BASIS_Q1_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int LL,LR,RR,RL,tot_parity,tot_ele,target_parity,count;
   int spin       = Model->spin_loc;
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int dim_onsite = Model->dim_onsite;
   int dim_LL     = System->Dim[LL_site];
   int dim_RR     = Enviro->Dim[RR_site];
   int *Dim       = GET_ARRAY_INT1(2);
   long dim_whole = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   //Calculate Dim
#pragma omp parallel for private (LR,RR,RL,tot_parity,tot_ele) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_parity = (System->Tot_Parity[LL_site][LL] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Parity[RR_site][RR] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(RL, spin))%2;
               tot_ele    = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(RL, spin);
               if (tot_ele == Model->tot_ele) {
#pragma omp critical
                  {
                     Dim[tot_parity]++;
                  }
               }
            }
         }
      }
   }
   
   DMRG_WHOLE_BASIS_Q1 *Basis = malloc(sizeof(DMRG_WHOLE_BASIS_Q1));
   Basis->max_dim      = FIND_MAX_INT1(Dim, 2);
   Basis->LL_LLLRRRRL  = GET_ARRAY_SINT2(2, Basis->max_dim);
   Basis->LR_LLLRRRRL  = GET_ARRAY_SINT2(2, Basis->max_dim);
   Basis->RL_LLLRRRRL  = GET_ARRAY_SINT2(2, Basis->max_dim);
   Basis->RR_LLLRRRRL  = GET_ARRAY_SINT2(2, Basis->max_dim);
   Basis->Inv_LLLRRRRL = Dmrg_Basis->Inv_LLLRRRRL;
   Basis->dim_RR       = dim_RR;
   Basis->dim_onsite   = dim_onsite;
   
   if (Basis->max_dim <= 0 || Dim[Model->tot_parity] != Dmrg_Basis->dim_LLLRRRRL) {
      printf("Error in GET_WHOLE_BASIS_Q1_SUPERBLOCK\n");
      printf("max_dim=%d\n",Basis->max_dim);
      printf("%d !=? %d\n", Dim[Model->tot_parity], Dmrg_Basis->dim_LLLRRRRL);
      exit(1);
   }
   
#pragma omp parallel for num_threads (Model->p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }
   
   Dim[0] = 0;
   Dim[1] = 0;
   
   target_parity = (Model->tot_parity + 1)%2;
#pragma omp parallel for private (LR,RR,RL,tot_parity,tot_ele,count,inv_sup) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_parity = (System->Tot_Parity[LL_site][LL] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Parity[RR_site][RR] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(RL, spin))%2;
               tot_ele    = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(RL, spin);
               if (tot_ele == Model->tot_ele && tot_parity == target_parity) {
#pragma omp critical
                  {
                     count = Dim[tot_parity];
                     Basis->LL_LLLRRRRL[tot_parity][count] = LL;
                     Basis->LR_LLLRRRRL[tot_parity][count] = LR;
                     Basis->RL_LLLRRRRL[tot_parity][count] = RL;
                     Basis->RR_LLLRRRRL[tot_parity][count] = RR;
                     inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
                     Basis->Inv_LLLRRRRL[inv_sup] = count;
                     Dim[tot_parity]++;
                  }
               }
            }
         }
      }
   }
   
   int basis;
   target_parity = Model->tot_parity;
   
#pragma omp parallel for private (LL,LR,RR,RL,inv_sup) num_threads (Model->p_threads)
   for (basis = 0; basis < Dmrg_Basis->dim_LLLRRRRL; basis++) {
      LL = Dmrg_Basis->LL_LLLRRRRL[basis];
      LR = Dmrg_Basis->LR_LLLRRRRL[basis];
      RL = Dmrg_Basis->RL_LLLRRRRL[basis];
      RR = Dmrg_Basis->RR_LLLRRRRL[basis];
      Basis->LL_LLLRRRRL[target_parity][basis] = LL;
      Basis->LR_LLLRRRRL[target_parity][basis] = LR;
      Basis->RL_LLLRRRRL[target_parity][basis] = RL;
      Basis->RR_LLLRRRRL[target_parity][basis] = RR;
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }
   
   Dim[target_parity] = Dmrg_Basis->dim_LLLRRRRL;
   
   Basis->Dim = Dim;
   
   return Basis;
   
}
