#include "dmrg.h"

void DMRG_EXPECTATION_ONSITE(CRS1 **M_LL, CRS1 *M_On, CRS1 **M_RR, double *Out, double *Vec, double *Temp_V, int p_threads, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int site;
   int LL_site      = Dmrg_Status->LL_site;
   int RR_site      = Dmrg_Status->RR_site;
   int dim_LLLRRRRL = Dmrg_Status->dim_LLLRRRRL;
   
   //LL_site
   for (site = 0; site <= LL_site; site++) {
      DMRG_V_M_LL_Q0(M_LL[site], Vec, Temp_V, dim_LLLRRRRL, p_threads, Dmrg_Basis);
      Out[site] = INNER_PRODUCT(Vec, Temp_V, dim_LLLRRRRL, p_threads);
   }
   
   //LR_site
   DMRG_V_M_LR_Q0(M_On, Vec, Temp_V, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   Out[LL_site + 1] = INNER_PRODUCT(Vec, Temp_V, dim_LLLRRRRL, p_threads);
   
   //RL_site
   DMRG_V_M_RL_Q0(M_On, Vec, Temp_V, dim_LLLRRRRL, p_threads, Dmrg_Basis);
   Out[LL_site + 2] = INNER_PRODUCT(Vec, Temp_V, dim_LLLRRRRL, p_threads);
   
   for (site = RR_site; site >= 0; site--) {
      DMRG_V_M_RR_Q0(M_RR[site], Vec, Temp_V, dim_LLLRRRRL, p_threads, Dmrg_Basis);
      Out[RR_site + LL_site + 3 - site] = INNER_PRODUCT(Vec, Temp_V, dim_LLLRRRRL, p_threads);
   }
   
   
}
