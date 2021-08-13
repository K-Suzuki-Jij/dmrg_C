//
//  GET_BASIS.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/25.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

long *GET_BASIS(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time, int dim) {
   
   Time->make_basis = omp_get_wtime();
   
   int site,tot_ele,tot_parity,site_basis,thread_num;
   int dim_onsite    = Model->dim_onsite;
   int spin_loc      = Model->spin_loc;
   int target_ele    = Model->tot_ele;
   int target_parity = Model->tot_parity;
   long tot_basis;
   long tot_dim = pow(Model->dim_onsite, Model->tot_site);
   long *Basis  = GET_ARRAY_LINT1(dim);
   int count = 0;
   
#pragma omp parallel for private (tot_ele,tot_parity,thread_num,site,site_basis) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele    = 0;
      tot_parity = 0;
      thread_num = omp_get_thread_num();
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  = EXACT_FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
         tot_parity += ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
      }
      if (tot_ele == target_ele && (tot_parity)%2 == target_parity) {
#pragma omp critical
         {
            Basis[count] = tot_basis;
            count++;
         }
      }
   }
   
   if (count != dim) {
      printf("Error in GET_BASIS\n");
      exit(1);
   }
   
   QUICK_SORT_LINT1(Basis, 0, dim);
   
   Time->make_basis = omp_get_wtime() - Time->make_basis;
   
   return Basis;
   
}
