#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "exact.h"
#include "SML.h"

void EXACT_V_M_Q1(CRS1 *M_On, int qn_out, double *Vec, int qn_in, double *Out_Vec, int dim_onsite, int site, int p_threads, EXACT_WHOLE_BASIS_Q1 *W_Basis) {
   
   if (qn_out < 0 || qn_in < 0) {
      printf("Error in EXACT_V_M_Q1\n");
      printf("qn_out=%d,qn_in=%d\n", qn_out, qn_in);
      exit(1);
   }
   
   int dim_out   = W_Basis->Dim[qn_out];
   int dim_in    = W_Basis->Dim[qn_in];
   long dim_site = (long)pow(dim_onsite, site);
   long whole_a_basis,i,j,inv,whole_target_basis;
   int local_basis;
   double val;
   
#pragma omp parallel for private (whole_target_basis,local_basis,val,j,whole_a_basis,inv) num_threads (p_threads)
   for (i = 0; i < dim_out; i++) {
      whole_target_basis = W_Basis->Basis[qn_out][i];
      local_basis = EXACT_FIND_SITE_STATE(whole_target_basis, site, dim_onsite);
      val = 0;
      for (j = M_On->Row[local_basis]; j < M_On->Row[local_basis + 1]; j++) {
         whole_a_basis = whole_target_basis - (local_basis - M_On->Col[j])*dim_site;
         inv = BINARY_SEARCH_LINT1(W_Basis->Basis[qn_in], 0, dim_in, whole_a_basis);
         if (inv >= 0) {
            val = val + Vec[inv]*M_On->Val[j];
         }
      }
      Out_Vec[i] = val;
   }
   
}
