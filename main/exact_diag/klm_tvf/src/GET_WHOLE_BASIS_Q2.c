//
//  GET_WHOLE_BASIS_Q2.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/30.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

EXACT_WHOLE_BASIS_Q2 *GET_WHOLE_BASIS_Q2(MODEL_1DKLM_TVF *Model, EXACT_TIME *Time) {
   
   double start = omp_get_wtime();
   
   int site,tot_ele,tot_parity,site_basis,thread_num,del_Ne;
   int dim_onsite    = Model->dim_onsite;
   int spin_loc      = Model->spin_loc;
   long tot_basis;
   long tot_dim      = pow(Model->dim_onsite, Model->tot_site);
   long ***Temp_Dim  = GET_ARRAY_LINT3(Model->p_threads, 3, 2);
   
   if (tot_dim <= 0 || pow(Model->dim_onsite, Model->tot_site) > ULONG_MAX) {
      printf("GET_WHOLE_BASIS_Q2\n");
      exit(1);
   }
   
#pragma omp parallel for private (tot_ele,tot_parity,thread_num,site,site_basis,del_Ne) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele    = 0;
      tot_parity = 0;
      thread_num = omp_get_thread_num();
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  = EXACT_FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
         tot_parity += ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
      }
      del_Ne     = Model->tot_ele - tot_ele;
      tot_parity = tot_parity%2;
      if (del_Ne == 0 || del_Ne == 1 || del_Ne == 2) {
         Temp_Dim[thread_num][del_Ne][tot_parity]++;
      }
   }
   
   for (del_Ne = 0; del_Ne <= 2; del_Ne++) {
      for (tot_parity = 0; tot_parity <= 1; tot_parity++) {
         for (thread_num = 1; thread_num < Model->p_threads; thread_num++) {
            Temp_Dim[0][del_Ne][tot_parity] = Temp_Dim[0][del_Ne][tot_parity] + Temp_Dim[thread_num][del_Ne][tot_parity];
         }
      }
   }
   
   for (del_Ne = 0; del_Ne <= 2; del_Ne++) {
      for (tot_parity = 0; tot_parity <= 1; tot_parity++) {
         if (Temp_Dim[0][del_Ne][tot_parity] <= 0 || Temp_Dim[0][del_Ne][tot_parity] > INT_MAX) {
            printf("Error in GET_WHOLE_BASIS_Q2\n");
            exit(1);
         }
      }
   }
   
   EXACT_WHOLE_BASIS_Q2 *W_Basis = malloc(sizeof(EXACT_WHOLE_BASIS_Q2));
   W_Basis->Dim     = GET_ARRAY_INT2(3, 2);
   W_Basis->max_dim = (int)FIND_MAX_LINT2(Temp_Dim[0], 3, 2);
   W_Basis->Basis   = GET_ARRAY_LINT3(3, 2, W_Basis->max_dim);
   
#pragma omp parallel for private (tot_ele,tot_parity,del_Ne,site,site_basis) num_threads (Model->p_threads)
   for (tot_basis = 0; tot_basis < tot_dim; tot_basis++) {
      tot_ele    = 0;
      tot_parity = 0;
      for (site = 0; site < Model->tot_site; site++) {
         site_basis  = EXACT_FIND_SITE_STATE(tot_basis, site, dim_onsite);
         tot_ele    += ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
         tot_parity += ONSITE_FIND_SITE_PARITY_SPIN_INV_BASIS_KLM(site_basis, spin_loc);
      }
      tot_parity = tot_parity%2;
      del_Ne     = Model->tot_ele - tot_ele;
      if (del_Ne == 0 || del_Ne == 1 || del_Ne == 2) {
#pragma omp critical
         {
            W_Basis->Basis[del_Ne][tot_parity][W_Basis->Dim[del_Ne][tot_parity]] = tot_basis;
            W_Basis->Dim[del_Ne][tot_parity]++;
         }
      }
   }
   
#pragma omp parallel for private (tot_parity) num_threads (Model->p_threads)
   for (del_Ne = 0; del_Ne <= 2; del_Ne++) {
      for (tot_parity = 0; tot_parity <= 1; tot_parity++) {
         if(W_Basis->Dim[del_Ne][tot_parity] != Temp_Dim[0][del_Ne][tot_parity]) {
            printf("Error in GET_WHOLE_BASIS_Q2\n");
            exit(1);
         }
         QUICK_SORT_LINT1(W_Basis->Basis[del_Ne][tot_parity], 0, W_Basis->Dim[del_Ne][tot_parity]);
      }
   }
   
   FREE_ARRAY_LINT3(Temp_Dim, Model->p_threads, 3);
   
   Time->make_basis = Time->make_basis + omp_get_wtime() - start;
   
   return W_Basis;
   
}

