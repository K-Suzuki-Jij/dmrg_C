//
//  EXPECTATION_INTERSITE_Q2.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/11.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"


void EXPECTATION_INTERSITE_Q2(CRS1 *M_O, CRS1 *M_R, int start, int end, double *Out, double *Vec, double **T_Vec1, double **T_Vec2, int dim_onsite, int p_threads, EXACT_WHOLE_BASIS_Q2 *W_Basis, int tot_sz) {
   
   int sz_in       = tot_sz;
   int sz_map_in   = (1 - SIGN(tot_sz))/2;
   int sz_out1     = tot_sz - 2;
   int sz_out2     = tot_sz + 2;
   int sz_map_out1 = (1 - SIGN(sz_out1))/2;
   int sz_map_out2 = (1 - SIGN(sz_out2))/2;
   int dim_out1    = W_Basis->Dim[sz_map_out1][abs(sz_out1)];
   int dim_out2    = W_Basis->Dim[sz_map_out2][abs(sz_out2)];

   EXACT_V_M_Q2(M_O, sz_map_out1, abs(sz_out1), Vec, sz_map_in, abs(sz_in), T_Vec1[0], "No", NULL, dim_onsite, start, p_threads, W_Basis);
   EXACT_V_M_Q2(M_O, sz_map_out2, abs(sz_out2), Vec, sz_map_in, abs(sz_in), T_Vec1[1], "No", NULL, dim_onsite, start, p_threads, W_Basis);
   
   int site,r;
   for (site = start; site < end; site++) {
      r = site - start;
      EXACT_V_M_Q2(M_R, sz_map_out1, abs(sz_out1), Vec, sz_map_in, abs(sz_in), T_Vec2[0], "No", NULL, dim_onsite, site, p_threads, W_Basis);
      EXACT_V_M_Q2(M_R, sz_map_out2, abs(sz_out2), Vec, sz_map_in, abs(sz_in), T_Vec2[1], "No", NULL, dim_onsite, site, p_threads, W_Basis);
      Out[r] = INNER_PRODUCT(T_Vec1[0], T_Vec2[0], dim_out1, p_threads) + INNER_PRODUCT(T_Vec1[1], T_Vec2[1], dim_out2, p_threads);
   }
   
}
