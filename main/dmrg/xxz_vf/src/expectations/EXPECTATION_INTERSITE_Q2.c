//
//  EXPECTATION_INTERSITE_Q2.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/04.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_INTERSITE_Q2(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_sz, double **Temp_V1, double **Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int site,r;
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int sz_in       = tot_sz;
   int sz_map_in   = (1 - SIGN(tot_sz))/2;
   int sz_out1     = tot_sz - 2;
   int sz_out2     = tot_sz + 2;
   int sz_map_out1 = (1 - SIGN(sz_out1))/2;
   int sz_map_out2 = (1 - SIGN(sz_out2))/2;
   int dim_in      = Dmrg_W_Basis->Dim[sz_map_in][abs(sz_in)];
   int dim_out1    = Dmrg_W_Basis->Dim[sz_map_out1][abs(sz_out1)];
   int dim_out2    = Dmrg_W_Basis->Dim[sz_map_out2][abs(sz_out2)];
   
   if (origin > LL_site + 1) {
      return;
   }
   
   if (origin == LL_site + 1) {
      DMRG_V_M_LR_Q2(M_On, sz_map_out1, abs(sz_out1), Vec, NULL, "No", Temp_V1[0], p_threads, Dmrg_W_Basis);
      DMRG_V_M_LR_Q2(M_On, sz_map_out2, abs(sz_out2), Vec, NULL, "No", Temp_V1[1], p_threads, Dmrg_W_Basis);
   }
   else {
      DMRG_V_M_LL_Q2(M_LL[origin], sz_map_out1, abs(sz_out1), Vec, Temp_V1[0], p_threads, Dmrg_W_Basis);
      DMRG_V_M_LL_Q2(M_LL[origin], sz_map_out2, abs(sz_out2), Vec, Temp_V1[1], p_threads, Dmrg_W_Basis);
   }
   
   //O*O to O*LL_site
   for (site = origin; site <= LL_site; site++) {
      r = site - origin;
      DMRG_V_M_LL_Q2(M_CF[r], sz_map_in, abs(sz_in), Vec, Temp_V2[0], p_threads, Dmrg_W_Basis);
      Out[r] = INNER_PRODUCT(Vec, Temp_V2[0], dim_in, p_threads);
   }
   
   //O*LR_site
   DMRG_V_M_LR_Q2(M_On, sz_map_out1, abs(sz_out1), Vec, NULL, "No", Temp_V2[0], p_threads, Dmrg_W_Basis);
   DMRG_V_M_LR_Q2(M_On, sz_map_out2, abs(sz_out2), Vec, NULL, "No", Temp_V2[1], p_threads, Dmrg_W_Basis);
   Out[LL_site + 1 - origin] = INNER_PRODUCT(Temp_V1[0], Temp_V2[0], dim_out1, p_threads);
   Out[LL_site + 1 - origin] = Out[LL_site + 1 - origin] + INNER_PRODUCT(Temp_V1[1], Temp_V2[1], dim_out2, p_threads);
   
   //O*RL_site
   DMRG_V_M_RL_Q2(M_On, sz_map_out1, abs(sz_out1), Vec, NULL, NULL, NULL, "No", Temp_V2[0], p_threads, Dmrg_W_Basis);
   DMRG_V_M_RL_Q2(M_On, sz_map_out2, abs(sz_out2), Vec, NULL, NULL, NULL, "No", Temp_V2[1], p_threads, Dmrg_W_Basis);
   Out[LL_site + 2 - origin] = INNER_PRODUCT(Temp_V1[0], Temp_V2[0], dim_out1, p_threads);
   Out[LL_site + 2 - origin] = Out[LL_site + 2 - origin] + INNER_PRODUCT(Temp_V1[1], Temp_V2[1], dim_out2, p_threads);
   
   //O*RR_site
   for (site = RR_site; site >= 0; site--) {
      r = RR_site + LL_site + 3 - site - origin;
      DMRG_V_M_RR_Q2(M_RR[site], sz_map_out1, abs(sz_out1), Vec, NULL, NULL, "No", Temp_V2[0], p_threads, Dmrg_W_Basis);
      DMRG_V_M_RR_Q2(M_RR[site], sz_map_out2, abs(sz_out2), Vec, NULL, NULL, "No", Temp_V2[1], p_threads, Dmrg_W_Basis);
      Out[r] = INNER_PRODUCT(Temp_V1[0], Temp_V2[0], dim_out1, p_threads);
      Out[r] = Out[r] + INNER_PRODUCT(Temp_V1[1], Temp_V2[1], dim_out2, p_threads);
   }
   
}
