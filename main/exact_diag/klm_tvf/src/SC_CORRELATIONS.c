//
//  SC_CORRELATIONS.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/30.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"
#define Max(i,j) (fabs(i) > fabs(j) ? fabs(i) : fabs(j))


void SC_CORRELATIONS(int sc_p, int site_start, double *GS_Vec, SC_MAT_1DKLM_TVF *Sc_Mat, EXACT_WHOLE_BASIS_Q2 *W_Basis, MODEL_1DKLM_TVF *Model) {
   
   if (site_start + 1 >= Model->tot_site || site_start <= 0) {
      return;
   }
   
   MAKE_SC_MAT_BASIS_1DKLM_TVF(sc_p, Sc_Mat, Model);
   
   HAM_BOX *Ham_Box     = GET_HAM_BOX(Model);
   double **Vec_CCSL_1  = GET_ARRAY_DOUBLE2(Sc_Mat->dim_ccsl_1   , W_Basis->max_dim);
   double **Vec_CCSL_2  = GET_ARRAY_DOUBLE2(Sc_Mat->dim_ccsl_2   , W_Basis->max_dim);
   double **Vec_CSL     = GET_ARRAY_DOUBLE2(Model->dim_csl_onsite, W_Basis->max_dim);
   double **Vec_CSL_CSL = GET_ARRAY_DOUBLE2(Sc_Mat->dim_csl_csl  , W_Basis->max_dim);
   double **Vec_Temp_1  = GET_ARRAY_DOUBLE2(Model->p_threads     , W_Basis->max_dim);
   double **Vec_Temp_2  = GET_ARRAY_DOUBLE2(Model->p_threads     , W_Basis->max_dim);

   int *N_Ele             = GET_ARRAY_INT1(Model->dim_onsite);
   int tot_p              = Model->tot_parity;
   int tot_site           = Model->tot_site;
   int dim_onsite         = Model->dim_onsite;
   int i,j,ind1,ind2,p1,p2,ii,jj,site,thread_num,r;
   
   for (i = 0; i < Model->dim_onsite; i++) {
      N_Ele[i] = ONSITE_FIND_SITE_ELE_SPIN_INV_BASIS_KLM(i, Model->spin_loc);
   }

   //Onsite_1 Reference
#pragma omp parallel for private (ind1,p1) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_ccsl_1; i++) {
      ind1 = Sc_Mat->CCSL_1_Num[i];
      p1   = Sc_Mat->CCSL_Parity[ind1];
      EXACT_V_M_Q2(Ham_Box->CCSL_On[ind1], 2, (p1 + tot_p)%2, GS_Vec, 0, tot_p, Vec_CCSL_1[i], "No", N_Ele, dim_onsite, site_start - 1, 1, W_Basis);
   }
   
   //Onsite_2 Reference
#pragma omp parallel for private (ind1,p1) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_ccsl_2; i++) {
      ind1 = Sc_Mat->CCSL_2_Num[i];
      p1   = Sc_Mat->CCSL_Parity[ind1];
      EXACT_V_M_Q2(Ham_Box->CCSL_On[ind1], 2, (p1 + tot_p)%2, GS_Vec, 0, tot_p, Vec_CCSL_2[i], "No", N_Ele, dim_onsite, site_start, 1, W_Basis);
   }

   //Intersite Reference 1
#pragma omp parallel for private (p1) num_threads (Model->p_threads)
   for (i = 0; i < Model->dim_csl_onsite; i++) {
      p1 = Sc_Mat->CSL_Parity[i];
      EXACT_V_M_Q2(Ham_Box->CSL_On[i], 1, (p1 + tot_p)%2, GS_Vec, 0, tot_p, Vec_CSL[i], "Yes", N_Ele, dim_onsite, site_start - 1, 1, W_Basis);
   }

   //Intersite Reference 1_2
#pragma omp parallel for private (ind1,ind2,p1,p2) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_csl_csl; i++) {
      ind1 = Sc_Mat->CSL_Num1[i];
      ind2 = Sc_Mat->CSL_Num2[i];
      p1   = Sc_Mat->CSL_Parity[ind1];
      p2   = Sc_Mat->CSL_Parity[ind2];
      EXACT_V_M_Q2(Ham_Box->CSL_On[ind1], 2, (p2 + p1 + tot_p)%2, Vec_CSL[ind2], 1, (p2 + tot_p)%2, Vec_CSL_CSL[i], "Yes", N_Ele, dim_onsite, site_start, 1, W_Basis);
   }

   //Reference * RR
#pragma omp parallel for private (thread_num,ind1,p1,site,r,j,jj,ii,ind2,p2) num_threads (Model->p_threads)
   for (i = 0; i < Sc_Mat->dim_tot; i++) {
      thread_num = omp_get_thread_num();
      
      //Length: Onsite_1
      if (i < Sc_Mat->dim_ccsl_1) {
         ind1 = Sc_Mat->CCSL_1_Num[i];
         p1   = Sc_Mat->CCSL_Parity[ind1];
         for (site = site_start; site < tot_site - 1; site++) {
            r = site - site_start;
            EXACT_V_M_Q2(Ham_Box->CCSL_On[ind1], 2, (p1 + tot_p)%2, GS_Vec, 0, tot_p, Vec_Temp_1[thread_num], "No", N_Ele, dim_onsite, site, 1, W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               
               //Reference: Onsite_2
               if (j < Sc_Mat->dim_ccsl_1) {
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_2[j], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_ccsl_1 <= j && j < Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl) {
                  jj = j - Sc_Mat->dim_ccsl_1;
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CSL_CSL[jj], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
               
               //Reference: Onsite_1
               else {
                 jj = j - (Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl);
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_1[jj], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
               
            }
         }
      }
      
      //Length: Intersite
      else if (Sc_Mat->dim_ccsl_1 <= i && i < Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl) {
         ii = i - Sc_Mat->dim_ccsl_1;
         ind1 = Sc_Mat->CSL_Num1[ii];
         ind2 = Sc_Mat->CSL_Num2[ii];
         p1   = Sc_Mat->CSL_Parity[ind1];
         p2   = Sc_Mat->CSL_Parity[ind2];
         for (site = site_start; site < tot_site - 1; site++) {
            r = site - site_start;
            EXACT_V_M_Q2(Ham_Box->CSL_On[ind2], 1, (p2 + tot_p)%2     , GS_Vec                , 0, tot_p         , Vec_Temp_1[thread_num], "Yes", N_Ele, dim_onsite, site + 1, 1, W_Basis);
            EXACT_V_M_Q2(Ham_Box->CSL_On[ind1], 2, (p1 + p2 + tot_p)%2, Vec_Temp_1[thread_num], 1, (p2 + tot_p)%2, Vec_Temp_2[thread_num], "Yes", N_Ele, dim_onsite, site    , 1, W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               
               //Reference: Onsite_2
               if (j < Sc_Mat->dim_ccsl_1) {
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_2[j], Vec_Temp_2[thread_num], W_Basis->Dim[2][(p1 + p2 + tot_p)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_ccsl_1 <= j && j < Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl) {
                  jj = j - Sc_Mat->dim_ccsl_1;
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CSL_CSL[jj], Vec_Temp_2[thread_num], W_Basis->Dim[2][(p1 + p2 + tot_p)%2], 1);
               }
               
               //Reference: Onsite_1
               else {
                  jj = j - (Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl);
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_1[jj], Vec_Temp_2[thread_num], W_Basis->Dim[2][(p1 + p2 + tot_p)%2], 1);
               }

            }
         }
      }
      
      //Length: Onsite_2
      else {
         ii   = i - (Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl);
         ind1 = Sc_Mat->CCSL_2_Num[ii];
         p1   = Sc_Mat->CCSL_Parity[ind1];
         for (site = site_start; site < tot_site - 1; site++) {
            r = site - site_start;
            EXACT_V_M_Q2(Ham_Box->CCSL_On[ind1], 2, (p1 + tot_p)%2, GS_Vec, 0, tot_p, Vec_Temp_1[thread_num], "No", N_Ele, dim_onsite, site + 1, 1, W_Basis);
            
            for (j = 0; j < Sc_Mat->dim_tot; j++) {
               
               //Reference: Onsite_2
               if (j < Sc_Mat->dim_ccsl_1) {
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_2[j], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
               //Reference: Intersite
               else if (Sc_Mat->dim_ccsl_1 <= j && j < Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl) {
                  jj = j - Sc_Mat->dim_ccsl_1;
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CSL_CSL[jj], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
               
               //Reference: Onsite_1
               else {
                  jj = j - (Sc_Mat->dim_ccsl_1 + Sc_Mat->dim_csl_csl);
                  Sc_Mat->Mat[r][j][i] = INNER_PRODUCT(Vec_CCSL_1[jj], Vec_Temp_1[thread_num], W_Basis->Dim[2][(p1 + tot_p)%2], 1);
               }
   
            }
         }
      }
   }
   
   FREE_ARRAY_DOUBLE2(Vec_CCSL_1   , Sc_Mat->dim_ccsl_1   );
   FREE_ARRAY_DOUBLE2(Vec_CCSL_2   , Sc_Mat->dim_ccsl_2   );
   FREE_ARRAY_DOUBLE2(Vec_CSL      , Model->dim_csl_onsite);
   FREE_ARRAY_DOUBLE2(Vec_CSL_CSL  , Sc_Mat->dim_csl_csl  );
   FREE_ARRAY_DOUBLE2(Vec_Temp_1   , Model->p_threads     );
   FREE_ARRAY_DOUBLE2(Vec_Temp_2   , Model->p_threads     );
   FREE_HAM_BOX(Ham_Box, Model);
   FREE_ARRAY_INT1(N_Ele);
   


}
