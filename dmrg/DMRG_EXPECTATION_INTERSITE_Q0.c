#include "dmrg.h"

void DMRG_EXPECTATION_INTERSITE_Q0(CRS1 **M_CF, CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, int origin, double *Vec, double *Temp_V1, double *Temp_V2, int p_threads, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int site,r;
   int LL_site      = Dmrg_Status->LL_site;
   int RR_site      = Dmrg_Status->RR_site;
   int dim_LLLRRRRL = Dmrg_Status->dim_LLLRRRRL;
   
   if (origin > LL_site + 1) {
      return;
   }
   
   if (origin == LL_site + 1) {
      DMRG_V_M_LR_Q0(M_On, Vec, Temp_V1, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   }
   else {
      DMRG_V_M_LL_Q0(M_LL[origin], Vec, Temp_V1, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   }
   
   //O*O to O*LL_site
   for (site = origin; site <= LL_site; site++) {
      r = site - origin;
      DMRG_V_M_LL_Q0(M_CF[r], Vec, Temp_V2, dim_LLLRRRRL, p_threads, Dmrg_Basis);
      Out[r] = INNER_PRODUCT(Vec, Temp_V2, dim_LLLRRRRL, p_threads);
   }
   
   //O*LR_site
   DMRG_V_M_LR_Q0(M_On, Vec, Temp_V2, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   Out[LL_site + 1 - origin] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_LLLRRRRL, p_threads);
   
   //O*RL_site
   DMRG_V_M_RL_Q0(M_On, Vec, Temp_V2, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   Out[LL_site + 2 - origin] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_LLLRRRRL, p_threads);
   
   //O*RR_site
   for (site = RR_site; site >= 0; site--) {
      r = RR_site + LL_site + 3 - site - origin;
      DMRG_V_M_RR_Q0(M_RR[site], Vec, Temp_V2, dim_LLLRRRRL, p_threads, Dmrg_Basis);
      Out[r] = INNER_PRODUCT(Temp_V1, Temp_V2, dim_LLLRRRRL, p_threads);
   }
   
}
