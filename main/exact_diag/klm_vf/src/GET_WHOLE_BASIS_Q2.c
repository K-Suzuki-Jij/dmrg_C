//
//  GET_WHOLE_BASIS_Q2.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/11.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

EXACT_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2(MODEL_1DKLM_VF *Model, EXACT_TIME *Time) {
   
   double start = omp_get_wtime();
   
   int site,tot_ele,tot_sz,site_basis,thread_num,sz_sign;
   int dim_onsite = Model->dim_onsite;
   int spin_loc   = Model->spin_loc;
   int target_ele = Model->tot_ele;
   int target_sz  = Model->tot_sz;
   int max_sz     = Model->tot_site*spin_loc + Model->tot_ele;
   long tot_basis;
   long tot_dim = pow(Model->dim_onsite, Model->tot_site);
   long ***Temp_Dim = GET_ARRAY_LINT3(Model->p_threads, 2, max_sz + 1);
   
   
   if (tot_dim <= 0 || pow(Model->dim_onsite, Model->tot_site) > ULONG_MAX) {
      printf("GET_WHOLE_BASIS_Q2 at 1\n");
      exit(1);
   }
   
#pragma omp parallel for private (tot_ele,tot_sz,thread_num,site,site_basis,sz_sign) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele    = 0;
      tot_sz = 0;
      thread_num = omp_get_thread_num();
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  = EXACT_FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SZBASIS_KLM(site_basis, spin_loc);
         tot_sz     += ONSITE_FIND_SITE_SZ_SZBASIS_KLM(site_basis, spin_loc);
      }
      sz_sign = SIGN(tot_sz);
      if (tot_ele == target_ele && (tot_sz == target_sz + 2 || tot_sz == target_sz || tot_sz == target_sz - 2)) {
         Temp_Dim[thread_num][(1 - sz_sign)/2][abs(tot_sz)]++;
      }
   }
   
   for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
      for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
         for (thread_num = 1; thread_num < Model->p_threads; thread_num++) {
            Temp_Dim[0][sz_sign][tot_sz] = Temp_Dim[0][sz_sign][tot_sz] + Temp_Dim[thread_num][sz_sign][tot_sz];
         }
      }
   }
   
   for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
      for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
         if (Temp_Dim[0][sz_sign][tot_sz] < 0 || Temp_Dim[0][sz_sign][tot_sz] > INT_MAX) {
            printf("Error in GET_WHOLE_BASIS_Q2 at 2\n");
            exit(1);
         }
      }
   }
   
   EXACT_WHOLE_BASIS_Q2 *W_Basis = malloc(sizeof(EXACT_WHOLE_BASIS_Q2));
   W_Basis->Dim     = GET_ARRAY_INT2(2, max_sz + 1);
   W_Basis->max_dim = (int)FIND_MAX_LINT2(Temp_Dim[0], 2, max_sz + 1);
   W_Basis->Basis   = malloc(sizeof(long**)*2);
   for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
      W_Basis->Basis[sz_sign] = malloc(sizeof(long*)*(max_sz + 1));
      for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
         W_Basis->Basis[sz_sign][tot_sz] = GET_ARRAY_LINT1(Temp_Dim[0][sz_sign][tot_sz]);
      }
   }
   
   
#pragma omp parallel for private (tot_ele,tot_sz,site,site_basis,sz_sign) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele = 0;
      tot_sz  = 0;
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  = EXACT_FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SZBASIS_KLM(site_basis, spin_loc);
         tot_sz     += ONSITE_FIND_SITE_SZ_SZBASIS_KLM(site_basis, spin_loc);
      }
      sz_sign = SIGN(tot_sz);
      
      if (tot_ele == target_ele && (tot_sz == target_sz + 2 || tot_sz == target_sz || tot_sz == target_sz - 2)) {
#pragma omp critical
         {
            W_Basis->Basis[(1 - sz_sign)/2][abs(tot_sz)][W_Basis->Dim[(1 - sz_sign)/2][abs(tot_sz)]] = tot_basis;
            W_Basis->Dim[(1 - sz_sign)/2][abs(tot_sz)]++;
         }
      }
   }
   

   
#pragma omp parallel for private (tot_sz) num_threads (Model->p_threads)
   for (sz_sign = 0; sz_sign <= 1; sz_sign++) {
      for (tot_sz = 0; tot_sz <= max_sz; tot_sz++) {
         if(W_Basis->Dim[sz_sign][tot_sz] != Temp_Dim[0][sz_sign][tot_sz]) {
            printf("Error in GET_WHOLE_BASIS_Q2 at 3\n");
            exit(1);
         }
         QUICK_SORT_LINT1(W_Basis->Basis[sz_sign][tot_sz], 0, W_Basis->Dim[sz_sign][tot_sz]);
      }
   }
   
   FREE_ARRAY_LINT3(Temp_Dim, Model->p_threads, 2);
   
   Time->make_basis = Time->make_basis + omp_get_wtime() - start;
   
   return W_Basis;
   
}
