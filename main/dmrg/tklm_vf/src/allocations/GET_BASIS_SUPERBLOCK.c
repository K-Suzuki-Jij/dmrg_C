//
//  GET_BASIS_SUPERBLOCK.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_BASIS *GET_BASIS_SUPERBLOCK(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_basis = omp_get_wtime();
   
   int LL,LR,RR,RL;
   int tot_sz,tot_ele_1,tot_ele_2,basis;
   int dim_LLLRRRRL = 0;
   
   int spin             = Model->spin_loc;
   int LL_site          = Dmrg_Status->LL_site;
   int RR_site          = Dmrg_Status->RR_site;
   int dim_onsite       = Model->dim_onsite;
   int dim_LL           = System->Dim[LL_site];
   int dim_RR           = Enviro->Dim[RR_site];
   int target_tot_ele_1 = (double)Model->tot_ele_1*(LL_site + RR_site + 4)/Model->tot_site;
   int target_tot_ele_2 = (double)Model->tot_ele_2*(LL_site + RR_site + 4)/Model->tot_site;
   int target_tot_sz    = (double)Model->tot_sz*(LL_site + RR_site + 4)/Model->tot_site;
   long dim_whole       = (long)dim_LL*dim_onsite*dim_RR*dim_onsite;
   long inv_sup;
   
   if (target_tot_ele_1 > target_tot_ele_2) {
      if ((target_tot_ele_1 + target_tot_ele_2)%2 == 1 && target_tot_sz%2 == 0) {
         target_tot_ele_1 = target_tot_ele_1 + 1;
      }
      if ((target_tot_ele_1 + target_tot_ele_2)%2 == 0 && target_tot_sz%2 == 1) {
         target_tot_ele_1 = target_tot_ele_1 + 1;
      }
   }
   else {
      if ((target_tot_ele_1 + target_tot_ele_2)%2 == 1 && target_tot_sz%2 == 0) {
         target_tot_ele_2 = target_tot_ele_2 + 1;
      }
      if ((target_tot_ele_1 + target_tot_ele_2)%2 == 0 && target_tot_sz%2 == 1) {
         target_tot_ele_2 = target_tot_ele_2 + 1;
      }
   }
      
   //Check Pint
   if (LL_site + RR_site + 4 == Model->tot_site) {
      if (target_tot_ele_1 != Model->tot_ele_1 || target_tot_ele_2 != Model->tot_ele_2 || target_tot_sz != Model->tot_sz) {
         printf("Error in GET_BASIS_SUPERBLOCK\n");
         printf("target_tot_ele_1=%d\n" ,target_tot_ele_1);
         printf("Model->tot_ele_1=%d\n" ,Model->tot_ele_1);
         printf("target_tot_ele_2=%d\n" ,target_tot_ele_2);
         printf("Model->tot_ele_2=%d\n" ,Model->tot_ele_2);
         printf("target_tot_sz=%1.1lf\n",(double)target_tot_sz/2);
         printf("Model->tot_sz=%1.1lf\n",(double)Model->tot_sz/2);
         exit(1);
      }
   }
   
   int temp_dim_LLLRRRRL = FIND_LLLRRRRL_DIM(target_tot_ele_1, target_tot_ele_2, target_tot_sz, System, Enviro, Model, Dmrg_Status);
      
   /* Allocate Basis */
   DMRG_BASIS *Basis = malloc(sizeof(DMRG_BASIS));
   DMRG_GET_BASIS_LLLRRRRL(temp_dim_LLLRRRRL, dim_onsite, dim_LL, dim_RR, Basis);
   
   short *Tot_Sz_LLLR_LLLRRRRL    = GET_ARRAY_SINT1(temp_dim_LLLRRRRL);
   short *Tot_Ele_1_LLLR_LLLRRRRL = GET_ARRAY_SINT1(temp_dim_LLLRRRRL);
   short *Tot_Ele_2_LLLR_LLLRRRRL = GET_ARRAY_SINT1(temp_dim_LLLRRRRL);
   
   dim_LLLRRRRL = 0;
#pragma omp parallel for private (LR,RR,RL,tot_sz,tot_ele_1,tot_ele_2) num_threads (Model->p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               tot_sz    = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Sz[RR_site][RR] + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(RL, spin);
               tot_ele_1 = System->Tot_Ele_1[LL_site][LL] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_1[RR_site][RR] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(RL, spin);
               tot_ele_2 = System->Tot_Ele_2[LL_site][LL] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(LR, spin) + Enviro->Tot_Ele_2[RR_site][RR] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(RL, spin);
               if (target_tot_sz == tot_sz && target_tot_ele_1 == tot_ele_1 && target_tot_ele_2 == tot_ele_2) {
#pragma omp critical
                  {
                     Basis->LL_LLLRRRRL[dim_LLLRRRRL]      = LL;
                     Basis->LR_LLLRRRRL[dim_LLLRRRRL]      = LR;
                     Basis->RL_LLLRRRRL[dim_LLLRRRRL]      = RL;
                     Basis->RR_LLLRRRRL[dim_LLLRRRRL]      = RR;
                     Tot_Sz_LLLR_LLLRRRRL[dim_LLLRRRRL]    = System->Tot_Sz[LL_site][LL] + ONSITE_FIND_SITE_SZ_SZBASIS_TKLM(LR, spin);
                     Tot_Ele_1_LLLR_LLLRRRRL[dim_LLLRRRRL] = System->Tot_Ele_1[LL_site][LL] + ONSITE_FIND_SITE_ELE_1_SZBASIS_TKLM(LR, spin);
                     Tot_Ele_2_LLLR_LLLRRRRL[dim_LLLRRRRL] = System->Tot_Ele_2[LL_site][LL] + ONSITE_FIND_SITE_ELE_2_SZBASIS_TKLM(LR, spin);
                     dim_LLLRRRRL++;
                  }
               }
            }
         }
      }
   }
   
   if (temp_dim_LLLRRRRL != dim_LLLRRRRL) {
      printf("Error in GET_BASIS_SUPERBLOCK\n");
      printf("dim_LLLRRRRL = %d, temp_dim_LLLRRRRL = %d\n", dim_LLLRRRRL, temp_dim_LLLRRRRL);
      exit(1);
   }
   
   DMRG_QUICK_SORT_BASIS_Q3(Basis->RL_LLLRRRRL, Basis->RR_LLLRRRRL, Basis->LR_LLLRRRRL, Basis->LL_LLLRRRRL, Tot_Sz_LLLR_LLLRRRRL, Tot_Ele_1_LLLR_LLLRRRRL, Tot_Ele_2_LLLR_LLLRRRRL, 0, dim_LLLRRRRL);
      
   Basis->dim_LLLRRRRL       = dim_LLLRRRRL;
   Basis->tot_sz_LLLRRRRL    = target_tot_sz;
   Basis->tot_ele_1_LLLRRRRL = target_tot_ele_1;
   Basis->tot_ele_2_LLLRRRRL = target_tot_ele_2;
   Basis->tot_ele_LLLRRRRL   = target_tot_ele_1 + target_tot_ele_2;
   Basis->dim_RR             = dim_RR;
   Basis->dim_LL             = dim_LL;
   Basis->dim_onsite         = dim_onsite;
   
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
   
   int dim_LLLR = 0;

   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      LL = Basis->LL_LLLRRRRL[basis];
      LR = Basis->LR_LLLRRRRL[basis];
      if (Basis->Inv_LLLR[LL][LR] == -1) {
         Basis->LL_LLLR[dim_LLLR]        = LL;
         Basis->LR_LLLR[dim_LLLR]        = LR;
         Basis->Tot_Sz_LLLR[dim_LLLR]    = Tot_Sz_LLLR_LLLRRRRL[basis];
         Basis->Tot_Ele_1_LLLR[dim_LLLR] = Tot_Ele_1_LLLR_LLLRRRRL[basis];
         Basis->Tot_Ele_2_LLLR[dim_LLLR] = Tot_Ele_2_LLLR_LLLRRRRL[basis];
         Basis->Tot_Ele_LLLR[dim_LLLR]   = Tot_Ele_1_LLLR_LLLRRRRL[basis] + Tot_Ele_2_LLLR_LLLRRRRL[basis];
         Basis->Inv_LLLR[LL][LR]         = dim_LLLR;
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
   Dmrg_Status->tot_ele_1    = target_tot_ele_1;
   Dmrg_Status->tot_ele_2    = target_tot_ele_2;
   Dmrg_Status->tot_ele      = target_tot_ele_1 + target_tot_ele_2;
   
   Dmrg_Status->dim_LL = System->Dim[Dmrg_Status->LL_site];
   Dmrg_Status->dim_RR = Enviro->Dim[Dmrg_Status->RR_site];
   
   FREE_ARRAY_SINT1(Tot_Sz_LLLR_LLLRRRRL);
   FREE_ARRAY_SINT1(Tot_Ele_1_LLLR_LLLRRRRL);
   FREE_ARRAY_SINT1(Tot_Ele_2_LLLR_LLLRRRRL);
   
   Dmrg_Time->make_basis = omp_get_wtime() - Dmrg_Time->make_basis;
      
   return Basis;
   
}
