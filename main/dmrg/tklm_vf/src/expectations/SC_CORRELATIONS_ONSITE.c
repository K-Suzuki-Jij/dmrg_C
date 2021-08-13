//
//  SC_CORRELATIONS_ONSITE.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/22.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void SC_CORRELATIONS_ONSITE(int sc_sz, int sc_ele_1, int sc_ele_2, SC_MAT_1DTKLM_VF *Sc_Mat, CRS1 **M_Reference, CRS1 ***M_Enviro, MODEL_1DTKLM_VF *Model, DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status, double *GS_Vec) {
   
   int RR_site    = Dmrg_Status->RR_site;
   int tot_sz     = Model->tot_sz;
   
   MAKE_SC_MAT_BASIS_1DTKLM_VF(sc_sz, sc_ele_1, sc_ele_2, Sc_Mat, Model);
   
   double **Vec_CCSL_RL = GET_ARRAY_DOUBLE2(Sc_Mat->dim_ccsl, Dmrg_W_Basis->max_dim);
   double **Vec_Temp    = GET_ARRAY_DOUBLE2(Model->p_threads, Dmrg_W_Basis->max_dim);
   
   int i,j,sz1,ind1,sz_map1,site,r,thread_num;
   
   //Onsite Reference RL
#pragma omp parallel for private (ind1,sz1,sz_map1) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_ccsl; i++) {
      ind1    = Sc_Mat->CCSL_Num[i];
      sz1     = Sc_Mat->CCSL_Sz[ind1] + tot_sz;
      sz_map1 = (1 - SIGN(sz1))/2;
      DMRG_V_M_RL_Q4(M_Reference[ind1], sc_ele_1, sc_ele_2, sz_map1, abs(sz1), GS_Vec, NULL, NULL, NULL, "No", Vec_CCSL_RL[i], Model->p_threads, Dmrg_W_Basis);
   }
   
#pragma omp parallel for private (thread_num,ind1,sz1,sz_map1,site,r,j) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_ccsl; i++) {
      thread_num  = omp_get_thread_num();
      ind1        = Sc_Mat->CCSL_Num[i];
      sz1         = Sc_Mat->CCSL_Sz[ind1] + tot_sz;
      sz_map1     = (1 - SIGN(sz1))/2;
      for (site = RR_site; site >= 0; site--) {
         r = RR_site - site + 1;
         for (j = 0; j < Sc_Mat->dim_ccsl; j++) {
            DMRG_V_M_RR_Q4(M_Enviro[ind1][site], sc_ele_1, sc_ele_2, sz_map1, abs(sz1), GS_Vec, NULL, NULL, "No", Vec_Temp[thread_num], 1, Dmrg_W_Basis);
            Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_RL[j], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[sc_ele_1][sc_ele_2][sz_map1][abs(sz1)], 1);
         }
      }
   }
   
   
   FREE_ARRAY_DOUBLE2(Vec_CCSL_RL, Sc_Mat->dim_ccsl);
   FREE_ARRAY_DOUBLE2(Vec_Temp   , Model->p_threads);
   
}
