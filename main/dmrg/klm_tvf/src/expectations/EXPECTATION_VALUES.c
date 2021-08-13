//
//  EXPECTATION_VALUES.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DKLM_TVF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status) {
   
   int c1 = (Dmrg_Status->sweep_now == 0);
   int c2 = (Dmrg_Status->sweep - Dmrg_Status->sweep_now >= 2);
   int c3 = (Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4 != Model->tot_site);
   int c4 = (Dmrg_Status->LL_site != Dmrg_Status->RR_site);
   
   if(c1 || c2 || c3 || c4){
      return;
   }
   
   OUTPUT_ENERGY(Model, Dmrg_Status);
   
   int tot_site     = Model->tot_site;
   int p_threads    = Model->p_threads;
   int cf_origin    = Model->cf_origin;
   int max_dim      = Dmrg_Param->max_dim_system;
   int elem_num     = max_dim*max_dim;
   int dim_LLLRRRRL = Dmrg_Status->dim_LLLRRRRL;
   
   CRS1 **M_Sys = GET_CRS2(Model->tot_site/2, max_dim, elem_num);
   CRS1 **M_Env = GET_CRS2(Model->tot_site/2, max_dim, elem_num);
   CRS1 **M_CF  = GET_CRS2(Model->tot_site, max_dim, elem_num);
   double *T_V1 = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   double *T_V2 = GET_ARRAY_DOUBLE1(dim_LLLRRRRL);
   
   ///SxL
   double *SxL = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SxL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SxL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SxL_On, M_Env, SxL, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SxL, "SxL", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxL, "Fourier_SxL", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SxL, "avg_SxL"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SxL_CF
   double *SxL_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->SxL_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->SxL_On, M_Env, SxL_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SxL_CF, SxL, "SxL_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxL_CF, "Fourier_SxL_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SxC
   double *SxC = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SxC_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SxC_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SxC_On, M_Env, SxC, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SxC, "SxC", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxC, "Fourier_SxC", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SxC, "avg_SxC"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SxC_CF
   double *SxC_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->SxC_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->SxC_On, M_Env, SxC_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SxC_CF, SxC, "SxC_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxC_CF, "Fourier_SxC_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///NC
   double *NC = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->NC_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->NC_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->NC_On, M_Env, NC, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(NC, "NC", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC, "Fourier_NC", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (NC, "avg_NC"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///NC_CF
   double *NC_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->NC_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->NC_On, M_Env, NC_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(NC_CF, NC, "NC_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC_CF, "Fourier_NC_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SCSL
   double *SCSL = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SCSL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SCSL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SCSL_On, M_Env, SCSL, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SCSL, "SCSL", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SCSL, "Fourier_SCSL", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SCSL, "avg_SCSL"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   FREE_ARRAY_DOUBLE1(T_V1);
   FREE_ARRAY_DOUBLE1(T_V2);
   
   DMRG_WHOLE_BASIS_Q1 *Dmrg_W_Basis = GET_WHOLE_BASIS_Q1_SUPERBLOCK(System, Enviro, Model, Dmrg_Basis, Dmrg_Status);
   
   double *TT_V1 = GET_ARRAY_DOUBLE1(Dmrg_W_Basis->max_dim);
   double *TT_V2 = GET_ARRAY_DOUBLE1(Dmrg_W_Basis->max_dim);
   
   ///SzL_CF
   double *SzL_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzL    = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SzL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SzL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_TWO(System->SzL_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   EXPECTATION_INTERSITE_PARITY(M_CF, M_Sys, System->SzL_On, M_Env, SzL_CF, Model->cf_origin, Vec, Model->tot_parity, TT_V1, TT_V2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SzL_CF, SzL, "SzL_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzL_CF, "Fourier_SzL_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SzC_CF
   double *SzC_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SzC    = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SzC_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SzC_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_TWO(System->SzC_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   EXPECTATION_INTERSITE_PARITY(M_CF, M_Sys, System->SzC_On, M_Env, SzC_CF, Model->cf_origin, Vec, Model->tot_parity, TT_V1, TT_V2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SzC_CF, SzC, "SzC_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzC_CF, "Fourier_SzC_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);

   FREE_ARRAY_DOUBLE1(TT_V1);
   FREE_ARRAY_DOUBLE1(TT_V2);
   
   FREE_WHOLE_BASIS_Q1_SUPERBLOCK(Dmrg_W_Basis);
   
   FREE_ARRAY_DOUBLE1(SzL);
   FREE_ARRAY_DOUBLE1(SzL_CF);
   FREE_ARRAY_DOUBLE1(SzC);
   FREE_ARRAY_DOUBLE1(SzC_CF);
   FREE_ARRAY_DOUBLE1(SxL);
   FREE_ARRAY_DOUBLE1(SxL_CF);
   FREE_ARRAY_DOUBLE1(SxC);
   FREE_ARRAY_DOUBLE1(SxC_CF);
   FREE_ARRAY_DOUBLE1(NC);
   FREE_ARRAY_DOUBLE1(NC_CF);
   FREE_ARRAY_DOUBLE1(SCSL);
   FREE_CRS2(M_Sys, Model->tot_site/2);
   FREE_CRS2(M_Env, Model->tot_site/2);
   FREE_CRS2(M_CF , Model->tot_site);
   
}
