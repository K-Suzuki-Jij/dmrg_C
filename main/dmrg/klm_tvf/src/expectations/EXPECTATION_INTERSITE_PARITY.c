//
//  EXPECTATION_INTERSITE_Q2.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/17.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_INTERSITE_PARITY(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, int tot_parity, double *Temp_V1, double *Temp_V2, int p_threads, DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int site,r;
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int parity_in  = tot_parity;
   int parity_out = (tot_parity + 1)%2;
   int dim_in     = Dmrg_W_Basis->Dim[parity_in];
   int dim_out    = Dmrg_W_Basis->Dim[parity_out];
   
   
   if (origin > LL_site + 1) {
      return;
   }
   
   if (origin == LL_site + 1) {
      DMRG_V_M_LR_Q1(M_On, parity_out, Vec, NULL, "No", Temp_V1, p_threads, Dmrg_W_Basis);
   }
   else {
      DMRG_V_M_LL_Q1(M_LL[origin], parity_out, Vec, Temp_V1, p_threads, Dmrg_W_Basis);
   }
   
   //O*O to O*LL_site
   for (site = origin; site <= LL_site; site++) {
      r = site - origin;
      DMRG_V_M_LL_Q1(M_CF[r], parity_in, Vec, Temp_V2, p_threads, Dmrg_W_Basis);
      Out[r] = INNER_PRODUCT(Vec, Temp_V2, dim_in, p_threads);
   }
   
   //O*LR_site
   DMRG_V_M_LR_Q1(M_On, parity_out, Vec, NULL, "No", Temp_V2, p_threads, Dmrg_W_Basis);
   Out[LL_site + 1 - origin] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_out, p_threads);
   
   //O*RL_site
   DMRG_V_M_RL_Q1(M_On, parity_out, Vec, NULL, NULL, NULL, "No", Temp_V2, p_threads, Dmrg_W_Basis);
   Out[LL_site + 2 - origin] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_out, p_threads);
   
   //O*RR_site
   for (site = RR_site; site >= 0; site--) {
      r = RR_site + LL_site + 3 - site - origin;
      DMRG_V_M_RR_Q1(M_RR[site], parity_out, Vec, NULL, NULL, "No", Temp_V2, p_threads, Dmrg_W_Basis);
      Out[r] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_out, p_threads);
   }
   
}
