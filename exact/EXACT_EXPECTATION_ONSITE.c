#include <stdio.h>
#include "exact.h"
#include "SML.h"

void EXACT_EXPECTATION_ONSITE(CRS1 *M_On, double *Out, double *Vec, double *T_Vec, int dim_onsite, int tot_site, int p_threads, EXACT_BASIS_INFO *Basis_Info) {
   
   int site;
   
   for (site = 0; site < tot_site; site++) {
      EXACT_V_M_Q0(M_On, Vec, T_Vec, dim_onsite, site, p_threads, Basis_Info);
      Out[site] = INNER_PRODUCT(Vec, T_Vec, Basis_Info->dim, p_threads);
   }
   
}

