#include <stdio.h>
#include "exact.h"
#include "SML.h"


void EXPECTATION_INTERSITE_Q1(CRS1 *M_O, CRS1 *M_R, int start, int end, double *Out, double *Vec, double *T_Vec1, double *T_Vec2, int dim_onsite, int p_threads, EXACT_WHOLE_BASIS_Q1 *W_Basis, int tot_parity) {
   
   int parity_in  = tot_parity;
   int parity_out = (tot_parity + 1)%2;
   int dim_out    = W_Basis->Dim[parity_out];
   
   EXACT_V_M_Q1(M_O, parity_out, Vec, parity_in, T_Vec1, dim_onsite, start, p_threads, W_Basis);
   
   int site,r;
   for (site = start; site < end; site++) {
      r = site - start;
      EXACT_V_M_Q1(M_R, parity_out, Vec, parity_in, T_Vec2, dim_onsite, site, p_threads, W_Basis);
      Out[r] = INNER_PRODUCT(T_Vec1, T_Vec2, dim_out, p_threads);
   }
   
}
