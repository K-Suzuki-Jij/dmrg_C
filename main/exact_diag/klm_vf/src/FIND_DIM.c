//
//  FIND_DIM.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int FIND_DIM(MODEL_1DKLM_VF *Model, EXACT_TIME *Time) {
   
   Time->find_dim = omp_get_wtime();
   
   int site,tot_ele,tot_sz,site_basis,thread_num;
   int dim_onsite = Model->dim_onsite;
   int spin_loc   = Model->spin_loc;
   int target_ele = Model->tot_ele;
   int target_sz  = Model->tot_sz;
   long tot_basis;
   long tot_dim   = pow(Model->dim_onsite, Model->tot_site);
   long *Temp_Dim = GET_ARRAY_LINT1(Model->p_threads);

   if (tot_dim <= 0 || pow(Model->dim_onsite, Model->tot_site) > ULONG_MAX) {
      printf("Error in FIND_DIM\n");
      exit(1);
   }
   
#pragma omp parallel for private (tot_ele,tot_sz,thread_num,site,site_basis) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele = 0;
      tot_sz  = 0;
      thread_num = omp_get_thread_num();
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  =  FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SZBASIS_KLM(site_basis, spin_loc);
         tot_sz     += ONSITE_FIND_SITE_SZ_SZBASIS_KLM(site_basis, spin_loc);
      }
      if (tot_ele == target_ele && tot_sz == target_sz) {
         Temp_Dim[thread_num]++;
      }
   }
   
   long dim = 0;
   
   for (thread_num = 0; thread_num < Model->p_threads; thread_num++) {
      dim = dim + Temp_Dim[thread_num];
   }
   
   if (dim <= 0 || dim > INT_MAX) {
      printf("Error in FIND_DIM\n");
      exit(1);
   }
   
   FREE_ARRAY_LINT1(Temp_Dim);
   
   Time->find_dim = omp_get_wtime() - Time->find_dim;
   
   return (int)dim;
   
}
