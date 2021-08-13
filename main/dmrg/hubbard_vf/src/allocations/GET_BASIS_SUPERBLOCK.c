//
//  GET_BASIS_SUPERBLOCK.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_basis = omp_get_wtime();
   
   int LL,LR,RR,RL;
   int tot_sz,tot_ele,basis;
   int dim_LLLRRRRL = 0;
   int dim_LLLR     = 0;
   
   int LL_site        = Dmrg_Status->LL_site;
   int RR_site        = Dmrg_Status->RR_site;
   int dim_onsite     = Model->dim_onsite;
   int dim_LL         = System->Dim[LL_site];
   int dim_RR         = Enviro->Dim[RR_site];
   int target_tot_ele = (double)Model->tot_ele*(LL_site + RR_site + 4)/Model->tot_site;
   int target_tot_sz  = (double)Model->tot_sz*(LL_site + RR_site + 4)/Model->tot_site;
   long dim_whole     = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   if (target_tot_ele%2 == 1 && target_tot_sz%2 == 0) {
      target_tot_ele = target_tot_ele + 1;
      //target_tot_sz  = target_tot_sz + 1;
   }
   if (target_tot_ele%2 == 0 && target_tot_sz%2 == 1) {
      target_tot_ele = target_tot_ele + 1;
      //target_tot_sz  = target_tot_sz + 1;
   }

   if (target_tot_ele == 0) {
      target_tot_ele = target_tot_ele + 1;
      target_tot_sz  = target_tot_sz + 1;
   }
   
   //Check Pint
   if (LL_site + RR_site + 4 == Model->tot_site) {
      if (target_tot_ele != Model->tot_ele || target_tot_sz != Model->tot_sz) {
         printf("Error in GET_BASIS_SUPERBLOCK\n");
         printf("target_tot_ele=%d\n",target_tot_ele);
         printf("Model->tot_el=%d\n" ,Model->tot_ele);
         printf("target_tot_sz=%1.1lf\n",(double)target_tot_sz/2);
         printf("Model->tot_sz=%1.1lf\n",(double)Model->tot_sz/2);
         exit(1);
      }
   }
   
   int temp_dim_LLLRRRRL = FIND_LLLRRRRL_DIM(target_tot_ele, target_tot_sz, System, Enviro, Model, Dmrg_Status);
   
   /* Allocate Basis */
   DMRG_BASIS *Basis = malloc(sizeof(DMRG_BASIS));
   DMRG_GET_BASIS_LLLRRRRL(temp_dim_LLLRRRRL, dim_onsite, dim_LL, dim_RR, Basis);
   
   short *Tot_Sz_LLLR_LLLRRRRL  = GET_ARRAY_SINT1(temp_dim_LLLRRRRL);
   short *Tot_Ele_LLLR_LLLRRRRL = GET_ARRAY_SINT1(temp_dim_LLLRRRRL);
   
   dim_LLLRRRRL = 0;
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz  = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_HUBBARD(LR) + Enviro->Tot_Sz[RR_site][RR] + ONSITE_FIND_SITE_SZ_SZBASIS_HUBBARD(RL);
               tot_ele = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD(LR) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD(RL);
               if (target_tot_sz == tot_sz && target_tot_ele == tot_ele) {
#pragma omp critical
                  {
                     Basis->LL_LLLRRRRL[dim_LLLRRRRL]    = LL;
                     Basis->LR_LLLRRRRL[dim_LLLRRRRL]    = LR;
                     Basis->RL_LLLRRRRL[dim_LLLRRRRL]    = RL;
                     Basis->RR_LLLRRRRL[dim_LLLRRRRL]    = RR;
                     Tot_Sz_LLLR_LLLRRRRL[dim_LLLRRRRL]  = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_HUBBARD(LR);
                     Tot_Ele_LLLR_LLLRRRRL[dim_LLLRRRRL] = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SZBASIS_HUBBARD(LR);
                     dim_LLLRRRRL++;
                  }
               }
            }
         }
      }
   }
   
   if (dim_LLLRRRRL <= 0) {
      printf("Error in GET_BASIS_SUPERBLOCK\n");
      printf("dim_LLLRRRRL = %d\n",dim_LLLRRRRL);
      exit(1);
   }
   
   DMRG_QUICK_SORT_BASIS_Q2(Basis->RL_LLLRRRRL, Basis->RR_LLLRRRRL, Basis->LR_LLLRRRRL, Basis->LL_LLLRRRRL, Tot_Sz_LLLR_LLLRRRRL, Tot_Ele_LLLR_LLLRRRRL, 0, dim_LLLRRRRL);
   
   Basis->dim_LLLRRRRL     = dim_LLLRRRRL;
   Basis->tot_sz_LLLRRRRL  = target_tot_sz;
   Basis->tot_ele_LLLRRRRL = target_tot_ele;
   Basis->dim_RR           = dim_RR;
   Basis->dim_LL           = dim_LL;
   Basis->dim_onsite       = dim_onsite;
   
#pragma omp parallel for num_threads (Model->p_threads)
   for (inv_sup = 0; inv_sup < dim_whole; inv_sup++) {
      Basis->Inv_LLLRRRRL[inv_sup] = -1;
   }

#pragma omp parallel for private (LL,LR,RL,RR,inv_sup) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      LL = Basis->LL_LLLRRRRL[basis];
      LR = Basis->LR_LLLRRRRL[basis];
      RL = Basis->RL_LLLRRRRL[basis];
      RR = Basis->RR_LLLRRRRL[basis];
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      Basis->Inv_LLLRRRRL[inv_sup] = basis;
   }
   
   DMRG_GET_BASIS_LLLR(dim_onsite, dim_LL, Basis);
   
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         Basis->Inv_LLLR[LL][LR] = -1;
      }
   }
   
   dim_LLLR = 0;
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      LL = Basis->LL_LLLRRRRL[basis];
      LR = Basis->LR_LLLRRRRL[basis];
      if (Basis->Inv_LLLR[LL][LR] == -1) {
         Basis->LL_LLLR[dim_LLLR]      = LL;
         Basis->LR_LLLR[dim_LLLR]      = LR;
         Basis->Tot_Sz_LLLR[dim_LLLR]  = Tot_Sz_LLLR_LLLRRRRL[basis];
         Basis->Tot_Ele_LLLR[dim_LLLR] = Tot_Ele_LLLR_LLLRRRRL[basis];
         Basis->Inv_LLLR[LL][LR]       = dim_LLLR;
         Basis->Sum_LLLR[dim_LLLR]++;
         dim_LLLR++;
      }
      else {
         Basis->Sum_LLLR[Basis->Inv_LLLR[LL][LR]]++;
      }
   }
   
   if (dim_LLLR <= 0) {
      printf("Error in GET_BASIS_SUPERBLOCK\n");
      printf("dim_LLLR = %d\n",dim_LLLR);
      exit(1);
   }
   
   Basis->dim_LLLR = dim_LLLR;
   
   DMRG_GET_BASIS_RRRL(dim_RR, dim_onsite, Basis);
   DMRG_GET_BASIS_LRRL(dim_onsite, Basis);
   
   Dmrg_Status->dim_LLLRRRRL = dim_LLLRRRRL;
   Dmrg_Status->tot_sz       = target_tot_sz;
   Dmrg_Status->tot_ele      = target_tot_ele;
   
   Dmrg_Status->dim_LL = System->Dim[Dmrg_Status->LL_site];
   Dmrg_Status->dim_RR = Enviro->Dim[Dmrg_Status->RR_site];
   
   FREE_ARRAY_SINT1(Tot_Sz_LLLR_LLLRRRRL);
   FREE_ARRAY_SINT1(Tot_Ele_LLLR_LLLRRRRL);
   
   Dmrg_Time->make_basis = omp_get_wtime() - Dmrg_Time->make_basis;
   
   return Basis;
   
}
