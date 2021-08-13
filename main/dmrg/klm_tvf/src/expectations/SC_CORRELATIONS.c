//
//  SC_CORRELATIONS.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/23.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void SC_CORRELATIONS(int sc_p, char LS_Couple[100], BLOCK *System, BLOCK *Enviro, SC_MAT_1DKLM_TVF *Sc_Mat, MODEL_1DKLM_TVF *Model, DMRG_WHOLE_BASIS_Q2 *Dmrg_W_Basis, DMRG_STATUS *Dmrg_Status, double *GS_Vec) {
   
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int tot_parity = Model->tot_parity;
   int dim_c_onsite;
   
   if (strcmp(LS_Couple, "Yes") == 0) {
      MAKE_LSC_MAT_BASIS_1DKLM_TVF(sc_p, Sc_Mat, Model->spin_loc, Model->tot_site);
      dim_c_onsite = (Model->spin_loc + 1)*(Model->spin_loc + 1)*2;
   }
   else {
      MAKE_SC_MAT_BASIS_1DKLM_TVF(sc_p, Sc_Mat, Model->tot_site);
      dim_c_onsite = 2;
   }
   
   double **Vec_CC_LR     = GET_ARRAY_DOUBLE2(Sc_Mat->dim_cc_1, Dmrg_W_Basis->max_dim);
   double **Vec_CC_RL     = GET_ARRAY_DOUBLE2(Sc_Mat->dim_cc_2, Dmrg_W_Basis->max_dim);
   double **Vec_C_LR      = GET_ARRAY_DOUBLE2(dim_c_onsite    , Dmrg_W_Basis->max_dim);
   double **Vec_C_C_LR_RL = GET_ARRAY_DOUBLE2(Sc_Mat->dim_c_c , Dmrg_W_Basis->max_dim);
   double **Vec_Temp      = GET_ARRAY_DOUBLE2(Model->p_threads, Dmrg_W_Basis->max_dim);
   int i,j,ind1,ind2,p1,p2,ii,jj,site,thread_num,r;
   
   //Onsite Reference LR
#pragma omp parallel for private (ind1,p1) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_cc_1; i++) {
      ind1 = Sc_Mat->CC_1_Num[i];
      p1   = (Sc_Mat->CC_Parity[ind1] + tot_parity)%2;
      DMRG_V_M_LR_Q2(Sc_Mat->CC_Onsite[ind1], 2, p1, GS_Vec, NULL, "No", Vec_CC_LR[i], 1, Dmrg_W_Basis);
   }
   
   //Onsite Reference RL
#pragma omp parallel for private (ind1,p1) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_cc_2; i++) {
      ind1 = Sc_Mat->CC_2_Num[i];
      p1   = (Sc_Mat->CC_Parity[ind1] + tot_parity)%2;
      DMRG_V_M_RL_Q2(Sc_Mat->CC_Onsite[ind1], 2, p1, GS_Vec, NULL, NULL, NULL, "No", Vec_CC_RL[i], 1, Dmrg_W_Basis);
   }
   
   //Intersite Reference LR
#pragma omp parallel for private (p1) num_threads (Model->p_threads)
   for (i = 0; i < dim_c_onsite; i++) {
      p1 = (Sc_Mat->C_Parity[i] + tot_parity)%2;
      DMRG_V_M_LR_Q2(Sc_Mat->C_Onsite[i], 1, p1, GS_Vec, System->Tot_Ele[LL_site], "Yes", Vec_C_LR[i], 1, Dmrg_W_Basis);
   }
   
   //Intersite Reference LR_RL
#pragma omp parallel for private (ind1,ind2,p1,p2) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_c_c; i++) {
      ind1 = Sc_Mat->C_Num1[i];
      ind2 = Sc_Mat->C_Num2[i];
      p1   = Sc_Mat->C_Parity[ind1];
      p2   = Sc_Mat->C_Parity[ind2];
      DMRG_V_M_RL_Q2(Sc_Mat->C_Onsite[ind1], 2, (p2 + p1 + tot_parity)%2, Vec_C_LR[ind2], System->Tot_Ele[LL_site], System->Tot_Ele[0], Enviro->Tot_Ele[RR_site], "Yes", Vec_C_C_LR_RL[i], 1, Dmrg_W_Basis);
   }
   
   //Reference * RR
#pragma omp parallel for private (thread_num,ind1,p1,site,r,j,jj,ii,ind2,p2) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_tot; i++) {
      thread_num = omp_get_thread_num();
      
      //Length: Onsite_r
      if (i < Sc_Mat->dim_cc_1) {
         ind1 = Sc_Mat->CC_1_Num[i];
         p1   = Sc_Mat->CC_Parity[ind1];
         for (site = RR_site; site >= 1; site--) {
            r = RR_site - site + 1;
            DMRG_V_M_RR_Q2(Sc_Mat->CC_Enviro[ind1][site], 2, (p1 + tot_parity)%2, GS_Vec, NULL, NULL, "No", Vec_Temp[thread_num], 1, Dmrg_W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               //Reference: Onsite_RL
               if (j < Sc_Mat->dim_cc_1) {
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CC_RL[j], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_cc_1 <= j && j < Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c) {
                  jj = j - Sc_Mat->dim_cc_1;
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_C_C_LR_RL[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
               //Reference: Onsite_LR
               else {
                  jj = j - (Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c);
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CC_LR[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
            }
         }
      }
      
      //Length: Intersite_r_r+1
      else if (Sc_Mat->dim_cc_1 <= i && i < Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c) {
         ii   = i - Sc_Mat->dim_cc_1;
         ind1 = Sc_Mat->C_Num1[ii];
         ind2 = Sc_Mat->C_Num2[ii];
         p1   = Sc_Mat->C_Parity[ind1];
         p2   = Sc_Mat->C_Parity[ind2];
         for (site = RR_site - 1; site >= 0; site--) {
            r = RR_site - 1 - site + 1;
            DMRG_V_M_RR_Q2(Sc_Mat->C_C_Enviro[ind2][ind1][site], 2, (p1 + p2 + tot_parity)%2, GS_Vec, NULL, NULL, "No", Vec_Temp[thread_num], 1, Dmrg_W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               //Reference: Onsite_RL
               if (j < Sc_Mat->dim_cc_1) {
                  Sc_Mat->Mat[r][j][i] = -1.0*INNER_PRODUCT(Vec_CC_RL[j], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + p2 + tot_parity)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_cc_1 <= j && j < Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c) {
                  jj = j - Sc_Mat->dim_cc_1;
                  Sc_Mat->Mat[r][j][i] = -1.0*INNER_PRODUCT(Vec_C_C_LR_RL[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + p2 + tot_parity)%2], 1);
               }
               //Reference: Onsite_LR
               else {
                  jj = j - (Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c);
                  Sc_Mat->Mat[r][j][i] = -1.0*INNER_PRODUCT(Vec_CC_LR[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + p2 + tot_parity)%2], 1);
               }
            }
         }
      }
      
      //Length: Onsite_r+1
      else {
         ii   = i - (Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c);
         ind1 = Sc_Mat->CC_2_Num[ii];
         p1   = Sc_Mat->CC_Parity[ind1];
         for (site = RR_site - 1; site >= 0; site--) {
            r = RR_site - 1 - site + 1;
            DMRG_V_M_RR_Q2(Sc_Mat->CC_Enviro[ind1][site], 2, (p1 + tot_parity)%2, GS_Vec, NULL, NULL, "No", Vec_Temp[thread_num], 1, Dmrg_W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               //Reference: Onsite_RL
               if (j < Sc_Mat->dim_cc_1) {
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CC_RL[j], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_cc_1 <= j && j < Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c) {
                  jj = j - Sc_Mat->dim_cc_1;
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_C_C_LR_RL[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
               //Reference: Onsite_LR
               else {
                  jj = j - (Sc_Mat->dim_cc_1 + Sc_Mat->dim_c_c);
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CC_LR[jj], Vec_Temp[thread_num], Dmrg_W_Basis->Dim[2][(p1 + tot_parity)%2], 1);
               }
            }
         }
      }
   }
   
   double f_norm;
   for (site = RR_site - 1; site >= 0; site--) {
      r = RR_site - 1 - site + 1;
      f_norm = 0;
      for (i = 0; i < Sc_Mat->dim_tot; i++) {
         for (j = 0; j < Sc_Mat->dim_tot; j++) {
            f_norm = f_norm + Sc_Mat->Mat[r][i][j]*Sc_Mat->Mat[r][i][j];
         }
      }
      Sc_Mat->F_Norm[r] = f_norm;
   }
   
   FREE_ARRAY_DOUBLE2(Vec_CC_LR    , Sc_Mat->dim_cc_1);
   FREE_ARRAY_DOUBLE2(Vec_CC_RL    , Sc_Mat->dim_cc_2);
   FREE_ARRAY_DOUBLE2(Vec_C_LR     , dim_c_onsite    );
   FREE_ARRAY_DOUBLE2(Vec_C_C_LR_RL, Sc_Mat->dim_c_c );
   FREE_ARRAY_DOUBLE2(Vec_Temp     , Model->p_threads);
   
}

