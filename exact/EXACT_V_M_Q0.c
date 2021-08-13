//
//  EXACT_V_M_Q0.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/26.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include <math.h>
#include <omp.h>
#include "exact.h"
#include "SML.h"

void EXACT_V_M_Q0(CRS1 *M_On, double *Vec, double *Out_Vec, int dim_onsite, int site, int p_threads, EXACT_BASIS_INFO *Basis_Info) {
   
   int dim       = Basis_Info->dim;
   long dim_site = (long)pow(dim_onsite, site);
   long whole_a_basis,i,j,inv,whole_target_basis;
   int local_basis;
   double val;
   
#pragma omp parallel for private (whole_target_basis,local_basis,val,j,whole_a_basis,inv) num_threads (p_threads)
   for (i = 0; i < dim; i++) {
      whole_target_basis = Basis_Info->Basis[i];
      local_basis = EXACT_FIND_SITE_STATE(whole_target_basis, site, dim_onsite);
      val = 0;
      for (j = M_On->Row[local_basis]; j < M_On->Row[local_basis + 1]; j++) {
         whole_a_basis = whole_target_basis - (local_basis - M_On->Col[j])*dim_site;
         inv = BINARY_SEARCH_LINT1(Basis_Info->Basis, 0, dim, whole_a_basis);
         if (inv >= 0) {
            val = val + Vec[inv]*M_On->Val[j];
         }
      }
      Out_Vec[i] = val;
   }

}
