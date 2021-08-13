//
//  EXPECTATION_SC_CORRELATIONS.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/17.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_SC_CORRELATIONS(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int c1 = (Dmrg_Status->sweep_now == 0);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c3 = (Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4 != Model->tot_site);
   int c4 = (Dmrg_Status->LL_site != Dmrg_Status->RR_site);
   int c5 = (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0);
   if(c1 || c2 || c3 || c4 || c5){
      return;
   }
   
   DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis = GET_WHOLE_BASIS_Q4_SUPERBLOCK(System, Enviro, Model, Dmrg_Basis, Dmrg_Status);
   SC_MAT_1DTKLM_VF    *Sc_Mat       = GET_SC_MAT_BASIS_1DTKLM_VF(Model);
   
   int del_sz;
   int max_del_sz = 4*Model->spin_loc + 2;
   int max_sz     = Model->spin_loc*Model->tot_site + Model->tot_ele_1 + Model->tot_ele_2;
   int min_sz     = -Model->spin_loc*Model->tot_site - (Model->tot_ele_1 + Model->tot_ele_2);
   char Name1[100],Name2[100];
   
   int tot_site        = Model->tot_site;
   int p_threads       = Model->p_threads;
   int spin_loc        = Model->spin_loc;
   int max_dim         = Dmrg_Param->max_dim_system;
   int dim_ccsl_onsite = Model->dim_ccsl_onsite;
   int dim_onsite      = Model->dim_onsite;
   int elem_num        = max_dim*max_dim*0.5;
   CRS1 **M_Reference  = GET_CRS2(dim_ccsl_onsite, dim_onsite, dim_onsite*dim_onsite);
   CRS1 ***M_Enviro    = GET_CRS3(dim_ccsl_onsite, Model->tot_site/2, max_dim, elem_num);
   
   int num;
   for (num = 0; num < dim_ccsl_onsite; num++) {
      ONSITE_CCSL_TKLM(num, spin_loc, M_Reference[num], 1);
      DMRG_TRANS_MAT_ONE(M_Reference[num], M_Enviro[num], Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   }
   
   for (del_sz = -max_del_sz; del_sz <= max_del_sz; del_sz++) {
      if (min_sz <= del_sz + Model->tot_sz && del_sz + Model->tot_sz <= max_sz) {
         SC_CORRELATIONS_ONSITE(del_sz, 2, 0, Sc_Mat, M_Reference, M_Enviro, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
         sprintf(Name1, "LSC_%d%d_%d_%d", 2, 0, del_sz, del_sz);
         sprintf(Name2, "%d%d_%d_%d", 2, 0, del_sz, del_sz);
         OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site, Name1, "ScCorrelations", Name2, Model, Dmrg_Status);
         
         SC_CORRELATIONS_ONSITE(del_sz, 1, 1, Sc_Mat, M_Reference, M_Enviro, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
         sprintf(Name1, "LSC_%d%d_%d_%d", 1, 1, del_sz, del_sz);
         sprintf(Name2, "%d%d_%d_%d", 1, 1, del_sz, del_sz);
         OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site, Name1, "ScCorrelations", Name2, Model, Dmrg_Status);
         
         SC_CORRELATIONS_ONSITE(del_sz, 0, 2, Sc_Mat, M_Reference, M_Enviro, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
         sprintf(Name1, "LSC_%d%d_%d_%d", 0, 2, del_sz, del_sz);
         sprintf(Name2, "%d%d_%d_%d", 0, 2, del_sz, del_sz);
         OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site, Name1, "ScCorrelations", Name2, Model, Dmrg_Status);
      }
   }
   
   FREE_CRS2(M_Reference, dim_ccsl_onsite);
   FREE_CRS3(M_Enviro, dim_ccsl_onsite, Model->tot_site/2);

   FREE_SC_MAT_BASIS_1DTKLM_VF(Sc_Mat, Model);
   FREE_WHOLE_BASIS_Q4_SUPERBLOCK(Dmrg_W_Basis, Model);
   
   FREE_T_MATRIX(System, Enviro, Model->tot_site, Dmrg_Status->sweep_now, Dmrg_Param->sweep, Dmrg_Param->Enviro_Copy);

}
