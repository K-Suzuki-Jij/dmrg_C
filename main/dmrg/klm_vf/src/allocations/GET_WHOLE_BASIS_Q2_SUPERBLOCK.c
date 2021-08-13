//
//  GET_WHOLE_BASIS_Q2_SUPERBLOCK.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/07.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int LL,LR,RR,RL,sz_sign,sz_map,tot_sz,tot_ele,count;
   int spin       = Model->spin_loc;
   int max_sz     = spin*Model->tot_site + Model->tot_ele + 4*spin + 2;
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int dim_onsite = Model->dim_onsite;
   int dim_LL     = System->Dim[LL_site];
   int dim_RR     = Enviro->Dim[RR_site];
   int **Dim      = GET_ARRAY_INT2(2, max_sz + 1);
   long dim_whole = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   //Calculate Dim
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele,sz_sign,sz_map) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz  = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_KLM(LR, spin) + Enviro->Tot_Sz[RR_site][RR] + ONSITE_FIND_SITE_SZ_SZBASIS_KLM(RL, spin);
               tot_ele = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SZBASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SZBASIS_KLM(RL, spin);
               if (tot_sz == Model->tot_sz || tot_sz == Model->tot_sz + 2 || tot_sz == Model->tot_sz - 2) {
                  if (tot_ele == Model->tot_ele) {
#pragma omp critical
                     {
                        sz_sign = SIGN(tot_sz);
                        sz_map  = (1 - sz_sign)/2;
                        Dim[sz_map][abs(tot_sz)]++;
                     }
                  }
               }
            }
         }
      }
   }
   
   
   DMRG_WHOLE_BASIS_Q2 *Basis = malloc(sizeof(DMRG_WHOLE_BASIS_Q2));
   Basis->max_dim      = FIND_MAX_INT2(Dim, 2, max_sz + 1);
   Basis->LL_LLLRRRRL  = malloc(sizeof(short**)*2);
   Basis->LR_LLLRRRRL  = malloc(sizeof(short**)*2);
   Basis->RL_LLLRRRRL  = malloc(sizeof(short**)*2);
   Basis->RR_LLLRRRRL  = malloc(sizeof(short**)*2);
   Basis->Inv_LLLRRRRL = Dmrg_Basis->Inv_LLLRRRRL;
   Basis->dim_RR       = dim_RR;
   Basis->dim_onsite   = dim_onsite;
   int i,j;
   for (i = 0; i < 2; i++) {
      Basis->LL_LLLRRRRL[i]  = malloc(sizeof(short*)*(max_sz + 1));
      Basis->LR_LLLRRRRL[i]  = malloc(sizeof(short*)*(max_sz + 1));
      Basis->RL_LLLRRRRL[i]  = malloc(sizeof(short*)*(max_sz + 1));
      Basis->RR_LLLRRRRL[i]  = malloc(sizeof(short*)*(max_sz + 1));
      for (j = 0; j < max_sz + 1; j++) {
         Basis->LL_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->LR_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->RL_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
         Basis->RR_LLLRRRRL[i][j]  = malloc(sizeof(short)*Dim[i][j]);
      }
   }
   
   if (Basis->max_dim <= 0) {
      printf("Error in GET_WHOLE_BASIS_Q2_SUPERBLOCK\n");
      printf("max_dim=%d\n",Basis->max_dim);
      exit(1);
   }
   
#pragma omp parallel for num_threads (Model->p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }
   
   for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
      for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
         Dim[sz_sign][tot_sz] = 0;
      }
   }
   
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele,sz_sign,sz_map,count,inv_sup) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_KLM(LR, spin) + Enviro->Tot_Sz[RR_site][RR] + ONSITE_FIND_SITE_SZ_SZBASIS_KLM(RL, spin);
               tot_ele = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SZBASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SZBASIS_KLM(RL, spin);
               if (tot_sz == Model->tot_sz || tot_sz == Model->tot_sz + 2 || tot_sz == Model->tot_sz - 2) {
                  if (tot_ele == Model->tot_ele) {
#pragma omp critical
                     {
                        sz_sign = SIGN(tot_sz);
                        sz_map  = (1 - sz_sign)/2;
                        count   = Dim[sz_map][abs(tot_sz)];
                        Basis->LL_LLLRRRRL[sz_map][abs(tot_sz)][count] = LL;
                        Basis->LR_LLLRRRRL[sz_map][abs(tot_sz)][count] = LR;
                        Basis->RL_LLLRRRRL[sz_map][abs(tot_sz)][count] = RL;
                        Basis->RR_LLLRRRRL[sz_map][abs(tot_sz)][count] = RR;
                        inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
                        Basis->Inv_LLLRRRRL[inv_sup] = count;
                        Dim[sz_map][abs(tot_sz)]++;
                     }
                  }
               }
            }
         }
      }
   }
   
   
   int basis;
   tot_sz = Dmrg_Basis->tot_sz_LLLRRRRL;
   sz_sign = SIGN(tot_sz);
   sz_map  = (1 - sz_sign)/2;
#pragma omp parallel for private (LL,LR,RR,RL,inv_sup) num_threads (Model->p_threads)
   for (basis = 0; basis < Dmrg_Basis->dim_LLLRRRRL; basis++) {
      LL     = Dmrg_Basis->LL_LLLRRRRL[basis];
      LR     = Dmrg_Basis->LR_LLLRRRRL[basis];
      RL     = Dmrg_Basis->RL_LLLRRRRL[basis];
      RR     = Dmrg_Basis->RR_LLLRRRRL[basis];
      Basis->LL_LLLRRRRL[sz_map][abs(tot_sz)][basis] = LL;
      Basis->LR_LLLRRRRL[sz_map][abs(tot_sz)][basis] = LR;
      Basis->RL_LLLRRRRL[sz_map][abs(tot_sz)][basis] = RL;
      Basis->RR_LLLRRRRL[sz_map][abs(tot_sz)][basis] = RR;
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }
   
   Basis->Dim = Dim;
   
   return Basis;
}
