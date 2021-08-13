#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "exact.h"
#include "onsite.h"
#include "SML.h"

void EXACT_V_M_Q3(CRS1 *M_On, int qn1_out, int qn2_out, int qn3_out, double *Vec, int qn1_in, int qn2_in, int qn3_in, double *Out_Vec, char Sign_Flag[], int *N_Ele, int dim_onsite, int op_site, int p_threads, EXACT_WHOLE_BASIS_Q3 *W_Basis) {
   
   if (qn1_out < 0 || qn2_out < 0 || qn3_out < 0 || qn1_in < 0 || qn2_in < 0 || qn3_in < 0) {
      printf("Error in EXACT_V_M_Q3\n");
      printf("qn1_out=%d,qn2_out=%d,qn3_out=%d,qn1_in=%d,qn2_in=%d,qn3_in=%d\n", qn1_out, qn2_out, qn3_out, qn1_in, qn2_in, qn3_in);
      exit(1);
   }
   
   int dim_out   = W_Basis->Dim[qn1_out][qn2_out][qn3_out];
   int dim_in    = W_Basis->Dim[qn1_in][qn2_in][qn3_in];
   long dim_site = (long)pow(dim_onsite, op_site);
   long whole_a_basis,i,j,inv,whole_target_basis;
   int local_basis,n_ele,site,sign;
   double val;
   
#pragma omp parallel for private (whole_target_basis,local_basis,n_ele,sign,val,j,whole_a_basis,inv) num_threads (p_threads)
   for (i = 0; i < dim_out; i++) {
      whole_target_basis = W_Basis->Basis[qn1_out][qn2_out][qn3_out][i];
      local_basis        = EXACT_FIND_SITE_STATE(whole_target_basis, op_site, dim_onsite);
      n_ele = 0;
      
      if (strcmp(Sign_Flag, "Yes") == 0) {
         for (site = 0; site < op_site; site++) {
            n_ele = n_ele + N_Ele[EXACT_FIND_SITE_STATE(whole_target_basis, site, dim_onsite)];
         }
         if (n_ele%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         sign = 1;
      }
      
      val = 0;
      for (j = M_On->Row[local_basis]; j < M_On->Row[local_basis + 1]; j++) {
         whole_a_basis = whole_target_basis - (local_basis - M_On->Col[j])*dim_site;
         inv = BINARY_SEARCH_LINT1(W_Basis->Basis[qn1_in][qn2_in][qn3_in], 0, dim_in, whole_a_basis);
         if (inv >= 0) {
            val = val + Vec[inv]*M_On->Val[j];
         }
         /*
         else {
            printf("Error in EXACT_V_M_Q2\n");
            exit(1);
         }
          */
      }
      Out_Vec[i] = val*sign;
   }
   
}
