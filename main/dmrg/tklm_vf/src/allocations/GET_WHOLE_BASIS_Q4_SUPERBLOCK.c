//
//  GET_WHOLE_BASIS_Q4_SUPERBLOCK.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/17.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_WHOLE_BASIS_Q4 *GET_WHOLE_BASIS_Q4_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int LL,LR,RR,RL,sz_sign,sz_map,tot_sz,tot_ele_1,tot_ele_2,count,del_Ne_1,del_Ne_2,c1,c2;
   int target_ele_1 = Model->tot_ele_1;
   int target_ele_2 = Model->tot_ele_2;
   int spin         = Model->spin_loc;
   int max_sz       = spin*Model->tot_site + Model->tot_ele_1 + Model->tot_ele_2 + 4*spin + 2;
   int LL_site      = Dmrg_Status->LL_site;
   int RR_site      = Dmrg_Status->RR_site;
   int dim_onsite   = Model->dim_onsite;
   int dim_LL       = System->Dim[LL_site];
   int dim_RR       = Enviro->Dim[RR_site];
   int ****Dim      = GET_ARRAY_INT4(3, 3, 2, max_sz + 1);
   long dim_whole   = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   //Calculate Dim
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele_1,tot_ele_2,del_Ne_1,del_Ne_2,c1,c2,sz_sign,sz_map) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz    = System->Tot_Sz[LL_site][LL]    + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(LR, spin)    + Enviro->Tot_Sz[RR_site][RR]    + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(RL, spin);
               tot_ele_1 = System->Tot_Ele_1[LL_site][LL] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_1[RR_site][RR] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(RL, spin);
               tot_ele_2 = System->Tot_Ele_2[LL_site][LL] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_2[RR_site][RR] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(RL, spin);
               del_Ne_1  = target_ele_1 - tot_ele_1;
               del_Ne_2  = target_ele_2 - tot_ele_2;
               c1 = (del_Ne_1 == 0 || del_Ne_1 == 1 || del_Ne_1 == 2);
               c2 = (del_Ne_2 == 0 || del_Ne_2 == 1 || del_Ne_2 == 2);
               
               if (c1 && c2) {
#pragma omp critical
                  {
                     sz_sign = SIGN(tot_sz);
                     sz_map  = (1 - sz_sign)/2;
                     Dim[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)]++;
                  }
               }
               
            }
         }
      }
   }
   
   DMRG_WHOLE_BASIS_Q4 *Basis = malloc(sizeof(DMRG_WHOLE_BASIS_Q4));
   Basis->max_dim      = FIND_MAX_INT4(Dim, 3, 3, 2, max_sz + 1);
   Basis->LL_LLLRRRRL  = malloc(sizeof(short****)*3);
   Basis->LR_LLLRRRRL  = malloc(sizeof(short****)*3);
   Basis->RL_LLLRRRRL  = malloc(sizeof(short****)*3);
   Basis->RR_LLLRRRRL  = malloc(sizeof(short****)*3);
   Basis->Inv_LLLRRRRL = Dmrg_Basis->Inv_LLLRRRRL;
   Basis->dim_RR       = dim_RR;
   Basis->dim_onsite   = dim_onsite;
   
   for (del_Ne_1 = 0; del_Ne_1 <= 2; del_Ne_1++) {
      Basis->LL_LLLRRRRL[del_Ne_1]  = malloc(sizeof(short***)*3);
      Basis->LR_LLLRRRRL[del_Ne_1]  = malloc(sizeof(short***)*3);
      Basis->RL_LLLRRRRL[del_Ne_1]  = malloc(sizeof(short***)*3);
      Basis->RR_LLLRRRRL[del_Ne_1]  = malloc(sizeof(short***)*3);
      
      for (del_Ne_2 = 0; del_Ne_2 <= 2; del_Ne_2++) {
         Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2]  = malloc(sizeof(short**)*2);
         Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2]  = malloc(sizeof(short**)*2);
         Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2]  = malloc(sizeof(short**)*2);
         Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2]  = malloc(sizeof(short**)*2);
         
         for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
            Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]  = malloc(sizeof(short*)*(max_sz + 1));
            Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]  = malloc(sizeof(short*)*(max_sz + 1));
            Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]  = malloc(sizeof(short*)*(max_sz + 1));
            Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign]  = malloc(sizeof(short*)*(max_sz + 1));
            
            for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
               Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz] = malloc(sizeof(short)*Dim[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz] = malloc(sizeof(short)*Dim[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz] = malloc(sizeof(short)*Dim[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
               Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_sign][tot_sz] = malloc(sizeof(short)*Dim[del_Ne_1][del_Ne_2][sz_sign][tot_sz]);
            }
         }
      }
   }
   
   if (Basis->max_dim <= 0) {
      printf("Error in GET_WHOLE_BASIS_Q4_SUPERBLOCK\n");
      printf("max_dim=%d\n",Basis->max_dim);
      exit(1);
   }
   
#pragma omp parallel for num_threads (Model->p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }
   
   
   for (del_Ne_1 = 0; del_Ne_1 <= 2; del_Ne_1++) {
      for (del_Ne_2 = 0; del_Ne_2 <= 2; del_Ne_2++) {
         for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
            for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
               Dim[del_Ne_1][del_Ne_2][sz_sign][tot_sz] = 0;
            }
         }
      }
   }
   
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele_1,tot_ele_2,del_Ne_1,del_Ne_2,c1,c2,sz_sign,sz_map,count,inv_sup) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz    = System->Tot_Sz[LL_site][LL]    + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(LR, spin)    + Enviro->Tot_Sz[RR_site][RR]    + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(RL, spin);
               tot_ele_1 = System->Tot_Ele_1[LL_site][LL] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_1[RR_site][RR] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(RL, spin);
               tot_ele_2 = System->Tot_Ele_2[LL_site][LL] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_2[RR_site][RR] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(RL, spin);
               del_Ne_1  = target_ele_1 - tot_ele_1;
               del_Ne_2  = target_ele_2 - tot_ele_2;
               c1 = (del_Ne_1 == 0 || del_Ne_1 == 1 || del_Ne_1 == 2);
               c2 = (del_Ne_2 == 0 || del_Ne_2 == 1 || del_Ne_2 == 2);
               
               if (c1 && c2) {
#pragma omp critical
                  {
                     sz_sign = SIGN(tot_sz);
                     sz_map  = (1 - sz_sign)/2;
                     count   = Dim[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)];
                     Basis->LL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)][count] = LL;
                     Basis->LR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)][count] = LR;
                     Basis->RL_LLLRRRRL[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)][count] = RL;
                     Basis->RR_LLLRRRRL[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)][count] = RR;
                     inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
                     Basis->Inv_LLLRRRRL[inv_sup] = count;
                     Dim[del_Ne_1][del_Ne_2][sz_map][abs(tot_sz)]++;
                  }
               }
            }
         }
      }
   }
   
   
   int basis;
   tot_sz  = Model->tot_sz;
   sz_sign = SIGN(tot_sz);
   sz_map  = (1 - sz_sign)/2;
#pragma omp parallel for private (LL,LR,RR,RL,inv_sup) num_threads (Model->p_threads)
   for (basis = 0; basis < Dmrg_Basis->dim_LLLRRRRL; basis++) {
      LL     = Dmrg_Basis->LL_LLLRRRRL[basis];
      LR     = Dmrg_Basis->LR_LLLRRRRL[basis];
      RL     = Dmrg_Basis->RL_LLLRRRRL[basis];
      RR     = Dmrg_Basis->RR_LLLRRRRL[basis];
      Basis->LL_LLLRRRRL[0][0][sz_map][abs(tot_sz)][basis] = LL;
      Basis->LR_LLLRRRRL[0][0][sz_map][abs(tot_sz)][basis] = LR;
      Basis->RL_LLLRRRRL[0][0][sz_map][abs(tot_sz)][basis] = RL;
      Basis->RR_LLLRRRRL[0][0][sz_map][abs(tot_sz)][basis] = RR;
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }
   
   Basis->Dim = Dim;
   
   return Basis;
}
