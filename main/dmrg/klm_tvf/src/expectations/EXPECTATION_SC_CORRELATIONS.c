//
//  EXPECTATION_SC_CORRELATIONS.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/22.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_SC_CORRELATIONS(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status) {
   
   int c1 = (Dmrg_Status->sweep_now == 0);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c3 = (Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4 != Model->tot_site);
   int c4 = (Dmrg_Status->LL_site != Dmrg_Status->RR_site);
   int c5 = (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0);
   if(c1 || c2 || c3 || c4 || c5){
      return;
   }
   
   DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis = GET_WHOLE_BASIS_Q2_SUPERBLOCK(System, Enviro, Model, Dmrg_Basis, Dmrg_Status);
   SC_MAT_1DKLM_TVF *Sc_Mat = GET_SC_MAT_BASIS_1DKLM_TVF(Model);
   
   int num,num1,num2;
   int tot_site   = Model->tot_site;
   int p_threads  = Model->p_threads;
   int spin_loc   = Model->spin_loc;
   int max_dim    = Dmrg_Param->max_dim_system;
   int elem_num   = max_dim*max_dim*0.3;
   int dim_onsite = Model->dim_onsite;
   int dim_lspin  = Model->dim_lspin;
   
   //Conventional SC
   int dim_cc_onsite  = 1;
   int dim_c_onsite   = 2;
   Sc_Mat->CC_Onsite  = GET_CRS2(dim_cc_onsite, dim_onsite, dim_onsite*dim_onsite);
   Sc_Mat->C_Onsite   = GET_CRS2(dim_c_onsite , dim_onsite, dim_onsite*dim_onsite);
   Sc_Mat->CC_Enviro  = GET_CRS3(dim_cc_onsite, Model->tot_site/2, max_dim, elem_num);
   Sc_Mat->C_C_Enviro = GET_CRS4(dim_c_onsite, dim_c_onsite, Model->tot_site/2, max_dim, elem_num);
   
   for (num1 = 0; num1 < dim_cc_onsite; num1++) {
      ONSITE_CC_KLM(num1, spin_loc, Sc_Mat->CC_Onsite[num1], 1.0);
      DMRG_TRANS_MAT_ONE(Sc_Mat->CC_Onsite[num1], Sc_Mat->CC_Enviro[num1], Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   }
   
   for (num1 = 0; num1 < dim_c_onsite; num1++) {
      ONSITE_C_KLM(num1, spin_loc, Sc_Mat->C_Onsite[num1], 1.0);
   }
   
#pragma omp parallel for private (num1,num2) num_threads (Model->p_threads)
   for (num = 0; num < dim_c_onsite*dim_c_onsite; num++) {
      num1 = num/dim_c_onsite;
      num2 = num%dim_c_onsite;
      DMRG_TRANS_MAT_C_C(Sc_Mat->C_Onsite[num1], Sc_Mat->C_Onsite[num2], Sc_Mat->C_C_Enviro[num1][num2], Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, Enviro->Tot_Ele, tot_site);
   }
   
   SC_CORRELATIONS(0, "No", System, Enviro, Sc_Mat, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site - 1, "SC_evev", "ScCorrelations", "EVEV", Model, Dmrg_Status);
   SC_CORRELATIONS(1, "No", System, Enviro, Sc_Mat, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site - 1, "SC_odod", "ScCorrelations", "ODOD", Model, Dmrg_Status);
   
   FREE_CRS2(Sc_Mat->CC_Onsite , dim_cc_onsite);
   FREE_CRS2(Sc_Mat->C_Onsite  , dim_c_onsite );
   FREE_CRS3(Sc_Mat->CC_Enviro , dim_cc_onsite, Model->tot_site/2);
   FREE_CRS4(Sc_Mat->C_C_Enviro, dim_c_onsite , dim_c_onsite, Model->tot_site/2);
   
   
   //Composite SC
   int dim_lscc_onsite = dim_lspin*dim_lspin;
   int dim_lsc_onsite  = dim_lspin*dim_lspin*2;
   Sc_Mat->CC_Onsite   = GET_CRS2(dim_lscc_onsite, dim_onsite, dim_onsite*dim_onsite);
   Sc_Mat->C_Onsite    = GET_CRS2(dim_lsc_onsite , dim_onsite, dim_onsite*dim_onsite);
   Sc_Mat->CC_Enviro   = GET_CRS3(dim_lscc_onsite, Model->tot_site/2, max_dim, elem_num);
   Sc_Mat->C_C_Enviro  = GET_CRS4(dim_lsc_onsite , dim_lsc_onsite, Model->tot_site/2, max_dim, elem_num);
   
   for (num1 = 0; num1 < dim_lscc_onsite; num1++) {
      ONSITE_LSCC_KLM(num1, spin_loc, Sc_Mat->CC_Onsite[num1], 1.0);
      DMRG_TRANS_MAT_ONE(Sc_Mat->CC_Onsite[num1], Sc_Mat->CC_Enviro[num1], Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   }
   
   for (num1 = 0; num1 < dim_lsc_onsite; num1++) {
      ONSITE_LSC_KLM(num1, spin_loc, Sc_Mat->C_Onsite[num1], 1.0);
   }
   
#pragma omp parallel for private (num1,num2) num_threads (Model->p_threads)
   for (num = 0; num < dim_lsc_onsite*dim_lsc_onsite; num++) {
      num1 = num/dim_lsc_onsite;
      num2 = num%dim_lsc_onsite;
      DMRG_TRANS_MAT_C_C(Sc_Mat->C_Onsite[num1], Sc_Mat->C_Onsite[num2], Sc_Mat->C_C_Enviro[num1][num2], Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, Enviro->Tot_Ele, tot_site);
   }
   
   SC_CORRELATIONS(0, "Yes", System, Enviro, Sc_Mat, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site - 1, "LSC_evev", "LScCorrelations", "EVEV", Model, Dmrg_Status);
   SC_CORRELATIONS(1, "Yes", System, Enviro, Sc_Mat, Model, Dmrg_W_Basis, Dmrg_Status, Vec);
   OUTPUT_SC_CORRELATIONS(Sc_Mat, Model->tot_site/2, Model->tot_site - 1, "LSC_odod", "LScCorrelations", "ODOD", Model, Dmrg_Status);
   
   FREE_CRS2(Sc_Mat->CC_Onsite , dim_lscc_onsite);
   FREE_CRS2(Sc_Mat->C_Onsite  , dim_lsc_onsite );
   FREE_CRS3(Sc_Mat->CC_Enviro , dim_lscc_onsite, Model->tot_site/2);
   FREE_CRS4(Sc_Mat->C_C_Enviro, dim_lsc_onsite , dim_lsc_onsite, Model->tot_site/2);
   
   FREE_WHOLE_BASIS_Q2_SUPERBLOCK(Dmrg_W_Basis);
   FREE_SC_MAT_BASIS_1DKLM_TVF(Sc_Mat, Model);
   
   FREE_T_MATRIX(System, Enviro, Model->tot_site, Dmrg_Status->sweep_now, Dmrg_Param->sweep, Dmrg_Param->Enviro_Copy);

}
