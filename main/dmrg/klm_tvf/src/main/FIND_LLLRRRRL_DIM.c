//
//  FIND_LLLRRRRL_DIM.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/11/05.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int FIND_LLLRRRRL_DIM(int target_tot_ele, int target_tot_parity, BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_STATUS *Dmrg_Status) {
   
   int LL,LR,RR,RL;
   int tot_parity,tot_ele,thread_num;
   int p_threads        = Model->p_threads;
   int spin             = Model->spin_loc;
   int LL_site          = Dmrg_Status->LL_site;
   int RR_site          = Dmrg_Status->RR_site;
   int dim_onsite       = Model->dim_onsite;
   int dim_LL           = System->Dim[LL_site];
   int dim_RR           = Enviro->Dim[RR_site];
   int *Dim_by_Threads  = GET_ARRAY_INT1(p_threads);
   
#pragma omp parallel for private (LR,RR,RL,thread_num,tot_parity,tot_ele) num_threads (p_threads)
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         for (RR = 0; RR < dim_RR; RR++) {
            for (RL = 0; RL < dim_onsite; RL++) {
               thread_num = omp_get_thread_num();
               tot_parity = (System->Tot_Parity[LL_site][LL] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Parity[RR_site][RR] + ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(RL, spin))%2;
               tot_ele    = System->Tot_Ele[LL_site][LL] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(LR, spin) + Enviro->Tot_Ele[RR_site][RR] + ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(RL, spin);
               if (target_tot_parity == tot_parity && target_tot_ele == tot_ele) {
                  Dim_by_Threads[thread_num]++;
               }
            }
         }
      }
   }
   
   int dim_LLLRRRRL = 0;
   
   for (thread_num = 0; thread_num < p_threads; thread_num++) {
      dim_LLLRRRRL = dim_LLLRRRRL + Dim_by_Threads[thread_num];
   }
   
   if (dim_LLLRRRRL <= 0 ) {
      printf("Error in FIND_LLLRRRRL_DIM\n");
      printf("dim_LLLRRRRL = %d\n",dim_LLLRRRRL);
      exit(1);
   }
   
   FREE_ARRAY_INT1(Dim_by_Threads);
   
   return dim_LLLRRRRL;
   
}
