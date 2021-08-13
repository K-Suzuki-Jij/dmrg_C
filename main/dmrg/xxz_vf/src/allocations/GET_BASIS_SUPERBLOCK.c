//
//  GET_BASIS_SUPERBLOCK.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_basis = omp_get_wtime();
   
   int LL,LR,RR,RL;
   int tot_sz,basis;
   int dim_LLLRRRRL = 0;
   
   int spin          = Model->spin;
   int target_tot_sz = Model->tot_sz;
   int LL_site       = Dmrg_Status->LL_site;
   int RR_site       = Dmrg_Status->RR_site;
   int dim_onsite    = Model->dim_onsite;
   int dim_LL        = System->Dim[LL_site];
   int dim_RR        = Enviro->Dim[RR_site];
   long dim_whole    = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   /* Allocate Basis */
   DMRG_BASIS *Basis = malloc(sizeof(DMRG_BASIS));
   DMRG_GET_BASIS_LLLRRRRL(dim_whole, dim_onsite, dim_LL, dim_RR, Basis);
   
   short *Tot_Sz_LLLR_LLLRRRRL = GET_ARRAY_SINT1(dim_whole);
   
   dim_LLLRRRRL = 0;
#pragma omp parallel for private (LR,RR,RL,tot_sz) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz = System->Tot_Sz[LL_site][LL] + FIND_SITE_SZ(LR, spin) + Enviro->Tot_Sz[RR_site][RR] + FIND_SITE_SZ(RL, spin);
               if (target_tot_sz == tot_sz) {
#pragma omp critical
                  {
                     Basis->LL_LLLRRRRL[dim_LLLRRRRL]   = LL;
                     Basis->LR_LLLRRRRL[dim_LLLRRRRL]   = LR;
                     Basis->RL_LLLRRRRL[dim_LLLRRRRL]   = RL;
                     Basis->RR_LLLRRRRL[dim_LLLRRRRL]   = RR;
                     Tot_Sz_LLLR_LLLRRRRL[dim_LLLRRRRL] = System->Tot_Sz[LL_site][LL] + FIND_SITE_SZ(LR, spin);
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
   
   DMRG_QUICK_SORT_BASIS_Q1(Basis->RL_LLLRRRRL, Basis->RR_LLLRRRRL, Basis->LR_LLLRRRRL, Basis->LL_LLLRRRRL, Tot_Sz_LLLR_LLLRRRRL, 0, dim_LLLRRRRL);
   
   Basis->dim_LLLRRRRL    = dim_LLLRRRRL;
   Basis->tot_sz_LLLRRRRL = target_tot_sz;
   Basis->dim_RR          = dim_RR;
   Basis->dim_LL          = dim_LL;
   Basis->dim_onsite      = dim_onsite;

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
   
   int dim_LLLR = 0;
   
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      LL = Basis->LL_LLLRRRRL[basis];
      LR = Basis->LR_LLLRRRRL[basis];
      if (Basis->Inv_LLLR[LL][LR] == -1) {
         Basis->LL_LLLR[dim_LLLR]     = LL;
         Basis->LR_LLLR[dim_LLLR]     = LR;
         Basis->Tot_Sz_LLLR[dim_LLLR] = Tot_Sz_LLLR_LLLRRRRL[basis];
         Basis->Inv_LLLR[LL][LR]      = dim_LLLR;
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
   
   Dmrg_Status->dim_LL = System->Dim[Dmrg_Status->LL_site];
   Dmrg_Status->dim_RR = Enviro->Dim[Dmrg_Status->RR_site];
   
   FREE_ARRAY_SINT1(Tot_Sz_LLLR_LLLRRRRL);
   
   Dmrg_Time->make_basis = omp_get_wtime() - Dmrg_Time->make_basis;
   
   return Basis;
   
}
