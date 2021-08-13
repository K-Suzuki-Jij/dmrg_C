//
//  EXPECTATION_VALUES.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void EXPECTATION_VALUES(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, double *Vec, DMRG_TIME *Dmrg_Time, DMRG_BASIS *Dmrg_Basis, DMRG_PARAMETER *Dmrg_Param, DMRG_STATUS *Dmrg_Status) {
   
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
   
   ///SzL
   double *SzL = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SzL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SzL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SzL_On, M_Env, SzL, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SzL, "SzL", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzL, "Fourier_SzL", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SzL, "avg_SzL"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SzL_CF
   double *SzL_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double s = omp_get_wtime();
   DMRG_TRANS_MAT_TWO(System->SzL_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   printf("SzLCF1=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->SzL_On, M_Env, SzL_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   printf("SzLCF2=%lf\n",omp_get_wtime() - s);
   OUTPUT_INTERSITE_VALUES(SzL_CF, SzL, "SzL_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzL_CF, "Fourier_SzL_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SzC_1
   double *SzC_1 = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SzC_1_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SzC_1_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SzC_1_On, M_Env, SzC_1, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SzC_1, "SzC_1", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzC_1, "Fourier_SzC_1", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SzC_1, "avg_SzC_1"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SzC_1_CF
   double *SzC_1_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->SzC_1_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->SzC_1_On, M_Env, SzC_1_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SzC_1_CF, SzC_1, "SzC_1_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzC_1_CF, "Fourier_SzC_1_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SzC_2
   double *SzC_2 = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SzC_2_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SzC_2_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SzC_2_On, M_Env, SzC_2, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SzC_2, "SzC_2", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzC_2, "Fourier_SzC_2", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SzC_2, "avg_SzC_2"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SzC_2_CF
   double *SzC_2_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->SzC_2_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->SzC_2_On, M_Env, SzC_2_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SzC_2_CF, SzC_2, "SzC_2_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SzC_2_CF, "Fourier_SzC_2_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);

   ///NC_1
   double *NC_1 = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->NC_1_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->NC_1_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->NC_1_On, M_Env, NC_1, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(NC_1, "NC_1", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC_1, "Fourier_NC_1", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (NC_1, "avg_NC_1"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///NC_1_CF
   double *NC_1_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->NC_1_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->NC_1_On, M_Env, NC_1_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(NC_1_CF, NC_1, "NC_1_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC_1_CF, "Fourier_NC_1_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///NC_2
   double *NC_2 = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->NC_2_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->NC_2_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->NC_2_On, M_Env, NC_2, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(NC_2, "NC_2", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC_2, "Fourier_NC_2", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (NC_2, "avg_NC_2"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///NC_2_CF
   double *NC_2_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_TWO(System->NC_2_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   DMRG_EXPECTATION_INTERSITE_Q0(M_CF, M_Sys, System->NC_2_On, M_Env, NC_2_CF, Model->cf_origin, Vec, T_V1, T_V2, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(NC_2_CF, NC_2, "NC_2_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(NC_2_CF, "Fourier_NC_2_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SC_1SL
   double *SC_1SL = GET_ARRAY_DOUBLE1(Model->tot_site);
   s = omp_get_wtime();
   DMRG_TRANS_MAT_ONE(System->SC_1SL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   printf("SCSL1=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_TRANS_MAT_ONE(Enviro->SC_1SL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   printf("SCSL2=%lf\n",omp_get_wtime() - s);
   s = omp_get_wtime();
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SC_1SL_On, M_Env, SC_1SL, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   printf("SCSL3=%lf\n",omp_get_wtime() - s);
   OUTPUT_ONSITE_VALUES(SC_1SL, "SC_1SL", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SC_1SL, "Fourier_SC_1SL", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SC_1SL, "avg_SC_1SL"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   ///SC_2SL
   double *SC_2SL = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SC_2SL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SC_2SL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_EXPECTATION_ONSITE(M_Sys, System->SC_2SL_On, M_Env, SC_2SL, Vec, T_V1, Model->p_threads, Dmrg_Basis, Dmrg_Status);
   OUTPUT_ONSITE_VALUES(SC_2SL, "SC_2SL", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SC_2SL, "Fourier_SC_2SL", 0, Model->tot_site, Model, Dmrg_Status);
   OUTPUT_AVERAGE_VALUES    (SC_2SL, "avg_SC_2SL"    , 0, Model->tot_site, Model, Dmrg_Status);
   
   FREE_ARRAY_DOUBLE1(T_V1);
   FREE_ARRAY_DOUBLE1(T_V2);
   
   DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis = GET_WHOLE_BASIS_Q2_SUPERBLOCK(System, Enviro, Model, Dmrg_Basis, Dmrg_Status);
   
   double **TT_V1 = GET_ARRAY_DOUBLE2(2, Dmrg_W_Basis->max_dim);
   double **TT_V2 = GET_ARRAY_DOUBLE2(2, Dmrg_W_Basis->max_dim);
   
   ///SxC_1_CF
   double *SxC_1_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxC_1    = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SxC_1_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SxC_1_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_TWO(System->SxC_1_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   EXPECTATION_INTERSITE_Q2(M_CF, M_Sys, System->SxC_1_On, M_Env, SxC_1_CF, Model->cf_origin, Vec, Model->tot_sz, TT_V1, TT_V2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SxC_1_CF, SxC_1, "SxC_1_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxC_1_CF, "Fourier_SxC_1_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SxC_2_CF
   double *SxC_2_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxC_2    = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SxC_2_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SxC_2_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_TWO(System->SxC_2_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   EXPECTATION_INTERSITE_Q2(M_CF, M_Sys, System->SxC_2_On, M_Env, SxC_2_CF, Model->cf_origin, Vec, Model->tot_sz, TT_V1, TT_V2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SxC_2_CF, SxC_2, "SxC_2_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxC_2_CF, "Fourier_SxC_2_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   ///SxL_CF
   double *SxL_CF = GET_ARRAY_DOUBLE1(Model->tot_site);
   double *SxL    = GET_ARRAY_DOUBLE1(Model->tot_site);
   DMRG_TRANS_MAT_ONE(System->SxL_On, M_Sys, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_ONE(Enviro->SxL_On, M_Env, Enviro->Dim, Enviro->Dim_LLLR, Enviro->TM, Enviro->TM_D, Enviro->Basis_LL_LLLR, Enviro->Basis_LR_LLLR, Enviro->Basis_Inv_LLLR, tot_site, p_threads);
   DMRG_TRANS_MAT_TWO(System->SxL_On, M_CF, System->Dim, System->Dim_LLLR, System->TM, System->TM_D, System->Basis_LL_LLLR, System->Basis_LR_LLLR, System->Basis_Inv_LLLR, cf_origin, tot_site, p_threads);
   EXPECTATION_INTERSITE_Q2(M_CF, M_Sys, System->SxL_On, M_Env, SxL_CF, Model->cf_origin, Vec, Model->tot_sz, TT_V1, TT_V2, Model->p_threads, Dmrg_W_Basis, Dmrg_Status);
   OUTPUT_INTERSITE_VALUES(SxL_CF, SxL, "SxL_CF", Model, Dmrg_Status);
   OUTPUT_FOURIER_COMPONENTS(SxL_CF, "Fourier_SxL_CF", 0, Model->tot_site - Model->cf_origin, Model, Dmrg_Status);
   
   FREE_ARRAY_DOUBLE2(TT_V1, 2);
   FREE_ARRAY_DOUBLE2(TT_V2, 2);
   
   FREE_WHOLE_BASIS_Q2_SUPERBLOCK(Dmrg_W_Basis, Model);
   
   FREE_ARRAY_DOUBLE1(SzL);
   FREE_ARRAY_DOUBLE1(SzL_CF);
   FREE_ARRAY_DOUBLE1(SzC_1);
   FREE_ARRAY_DOUBLE1(SzC_1_CF);
   FREE_ARRAY_DOUBLE1(SzC_2);
   FREE_ARRAY_DOUBLE1(SzC_2_CF);
   FREE_ARRAY_DOUBLE1(SxC_1);
   FREE_ARRAY_DOUBLE1(SxC_1_CF);
   FREE_ARRAY_DOUBLE1(SxC_2);
   FREE_ARRAY_DOUBLE1(SxC_2_CF);
   FREE_ARRAY_DOUBLE1(SxL);
   FREE_ARRAY_DOUBLE1(SxL_CF);
   FREE_ARRAY_DOUBLE1(NC_1);
   FREE_ARRAY_DOUBLE1(NC_1_CF);
   FREE_ARRAY_DOUBLE1(NC_2);
   FREE_ARRAY_DOUBLE1(NC_2_CF);
   FREE_ARRAY_DOUBLE1(SC_1SL);
   FREE_ARRAY_DOUBLE1(SC_2SL);
   FREE_CRS2(M_Sys, Model->tot_site/2);
   FREE_CRS2(M_Env, Model->tot_site/2);
   FREE_CRS2(M_CF , Model->tot_site);
   
}
