//
//  GET_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/23.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int LL,LR,RR,RL,tot_parity,tot_ele,target_parity,count,del_Ne;
   int spin       = Model->spin_loc;
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int dim_onsite = Model->dim_onsite;
   int dim_LL     = System->Dim[LL_site];
   int dim_RR     = Enviro->Dim[RR_site];
   int **Dim      = GET_ARRAY_INT2(3, 2);
   long dim_whole = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   //Calculate Dim
#pragma omp parallel for private (LR,RR,RL,tot_parity,tot_ele,del_Ne) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_parity = (System->Tot_Parity[LL_site][LL] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Parity[RR_site][RR] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(RL, spin))%2;
               tot_ele    = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(RL, spin);
               del_Ne     = Model->tot_ele - tot_ele;
               if (del_Ne == 0 || del_Ne == 1 || del_Ne == 2) {
#pragma omp critical
                  {
                     Dim[del_Ne][tot_parity]++;
                  }
               }
            }
         }
      }
   }
   
   DMRG_WHOLE_BASIS_Q2 *Basis = malloc(sizeof(DMRG_WHOLE_BASIS_Q2));
   Basis->max_dim      = FIND_MAX_INT2(Dim, 3, 2);
   Basis->LL_LLLRRRRL  = malloc(sizeof(short**)*3);
   Basis->LR_LLLRRRRL  = malloc(sizeof(short**)*3);
   Basis->RL_LLLRRRRL  = malloc(sizeof(short**)*3);
   Basis->RR_LLLRRRRL  = malloc(sizeof(short**)*3);
   Basis->Inv_LLLRRRRL = Dmrg_Basis->Inv_LLLRRRRL;
   Basis->dim_RR       = dim_RR;
   Basis->dim_onsite   = dim_onsite;
   int i,j;
   for (i = 0; i <= 2; i++) {
      Basis->LL_LLLRRRRL[i]  = malloc(sizeof(short*)*2);
      Basis->LR_LLLRRRRL[i]  = malloc(sizeof(short*)*2);
      Basis->RL_LLLRRRRL[i]  = malloc(sizeof(short*)*2);
      Basis->RR_LLLRRRRL[i]  = malloc(sizeof(short*)*2);
      for (j = 0; j <= 1; j++) {
         Basis->LL_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->LR_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->RL_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->RR_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
      }
   }
   
   
   if (Basis->max_dim <= 0 || Dim[0][Model->tot_parity] != Dmrg_Basis->dim_LLLRRRRL) {
      printf("Error in GET_WHOLE_BASIS_Q2_SUPERBLOCK\n");
      printf("max_dim=%d\n",Basis->max_dim);
      printf("%d !=? %d\n", Dim[0][Model->tot_parity], Dmrg_Basis->dim_LLLRRRRL);
      exit(1);
   }
   
#pragma omp parallel for num_threads (Model->p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }
   
   for (i = 0; i <= 2; i++) {
      for (j = 0; j <= 1; j++) {
         Dim[i][j] = 0;
      }
   }
   
#pragma omp parallel for private (LR,RR,RL,tot_parity,tot_ele,del_Ne,count,inv_sup) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_parity = (System->Tot_Parity[LL_site][LL] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Parity[RR_site][RR] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(RL, spin))%2;
               tot_ele    = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(RL, spin);
               del_Ne     = Model->tot_ele - tot_ele;
               if (del_Ne == 0 || del_Ne == 1 || del_Ne == 2) {
#pragma omp critical
                  {
                     count = Dim[del_Ne][tot_parity];
                     Basis->LL_LLLRRRRL[del_Ne][tot_parity][count] = LL;
                     Basis->LR_LLLRRRRL[del_Ne][tot_parity][count] = LR;
                     Basis->RL_LLLRRRRL[del_Ne][tot_parity][count] = RL;
                     Basis->RR_LLLRRRRL[del_Ne][tot_parity][count] = RR;
                     inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
                     Basis->Inv_LLLRRRRL[inv_sup] = count;
                     Dim[del_Ne][tot_parity]++;
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
      Basis->LL_LLLRRRRL[0][target_parity][basis] = LL;
      Basis->LR_LLLRRRRL[0][target_parity][basis] = LR;
      Basis->RL_LLLRRRRL[0][target_parity][basis] = RL;
      Basis->RR_LLLRRRRL[0][target_parity][basis] = RR;
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }
   
   Basis->Dim = Dim;
   
   return Basis;
   
}


