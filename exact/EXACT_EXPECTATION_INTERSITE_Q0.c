#include <stdio.h>
#include "exact.h"
#include "SML.h"


void EXACT_EXPECTATION_INTERSITE_Q0(CRS1 *M_O, CRS1 *M_R, int start, int end, double *Out, double *Vec, double *T_Vec1, double *T_Vec2, int dim_onsite, int p_threads, EXACT_BASIS_INFO *Basis_Info) {
   
   
   EXACT_V_M_Q0(M_O, Vec, T_Vec1, dim_onsite, start, p_threads, Basis_Info);
   
   int site,r;
   for (site = start; site < end; site++) {
      r = site - start;
      EXACT_V_M_Q0(M_R, Vec, T_Vec2, dim_onsite, site, p_threads, Basis_Info);
      Out[r] = INNER_PRODUCT(T_Vec1, T_Vec2, Basis_Info->dim, p_threads);
   }
   
}



