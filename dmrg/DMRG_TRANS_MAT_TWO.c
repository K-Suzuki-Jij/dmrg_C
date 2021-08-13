#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SML.h"
#include "dmrg.h"

void DMRG_TRANS_MAT_TWO(CRS1 *M_On, CRS1 **Out, int *Dim_LL, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int cf_origin, int tot_site, int p_threads) {
   
   int LL_site,r,thread_num;
   int max_elem_num_LLLR = 0;
   int max_elem_num_LL   = 0;
   int max_dim_LLLR      = 0;
   int max_dim_LL        = 0;

   for (LL_site = 0; LL_site <= tot_site/2 - 2; LL_site++) {
      if (max_elem_num_LLLR < Dim_LLLR[LL_site]*Dim_LL[LL_site + 1]) {
         max_elem_num_LLLR = Dim_LLLR[LL_site]*Dim_LL[LL_site + 1];
      }
      if (max_dim_LLLR < Dim_LLLR[LL_site]) {
         max_dim_LLLR = Dim_LLLR[LL_site];
      }
   }
   for (LL_site = 0; LL_site < tot_site; LL_site++) {
      if (max_dim_LL < Dim_LL[LL_site]) {
         max_dim_LL = Dim_LL[LL_site];
      }
   }
      
   max_elem_num_LL = max_dim_LL*max_dim_LL;
   
   CRS1 **M_CRS = GET_CRS2(p_threads, max_dim_LLLR, max_elem_num_LLLR);
   CCS1 **M_CCS = GET_CCS2(p_threads, max_dim_LLLR, max_elem_num_LLLR);
   CRS1 *M_Base = GET_CRS1(max_dim_LL, max_elem_num_LL);
   CRS1 *MM_On_T= GET_CRS1(max_dim_LL, max_elem_num_LL);
   CRS1 *MM_On  = GET_CRS1(max_dim_LL, max_elem_num_LL);
   double **V   = GET_ARRAY_DOUBLE2(p_threads, max_dim_LLLR);
   
   MATRIX_PRODUCT_CRS1(M_On, M_On, MM_On);
   
   DMRG_BASIS_LLLR *Basis = malloc(sizeof(DMRG_BASIS_LLLR));
   
   for (LL_site = 0; LL_site <= tot_site/2 - 3; LL_site++) {
      
      CCS1 *T_M       = TM[LL_site];
      CRS1 *T_MD      = TM_D[LL_site];
      Basis->dim_LLLR = Dim_LLLR[LL_site];
      Basis->LL_LLLR  = LL_LLLR[LL_site];
      Basis->LR_LLLR  = LR_LLLR[LL_site];
      Basis->Inv_LLLR = Inv_LLLR[LL_site];
      
      //Make Base Matrix
      if (cf_origin == 0 && LL_site == 0) {
         COPY_CRS1(M_On, M_Base, 1);
         COPY_CRS1(MM_On, MM_On_T, 1);
      }
      else if (cf_origin == LL_site + 1) {
         DMRG_EXTEND_AND_REDUCE(NULL, M_On, M_Base, NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
         DMRG_EXTEND_AND_REDUCE(NULL, MM_On, MM_On_T, NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
      }
            
      if (LL_site >= cf_origin) {
         r = LL_site - cf_origin;
         DMRG_EXTEND_AND_REDUCE(M_Base, M_On, Out[r + 1], NULL, "LL_LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
         
         if (r == 0) {
            DMRG_EXTEND_AND_REDUCE(MM_On_T, NULL, Out[r], NULL, "LL", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
         }
         else {
#pragma omp parallel for private (thread_num) num_threads (p_threads)
            for (r = 0; r <= LL_site - cf_origin; r++){
               thread_num = omp_get_thread_num();
               DMRG_EXTEND_AND_REDUCE(Out[r], NULL, Out[r], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Basis);
            }
         }
         
         DMRG_EXTEND_AND_REDUCE(M_Base , NULL, M_Base , NULL, "LL", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
         DMRG_EXTEND_AND_REDUCE(MM_On_T, NULL, MM_On_T, NULL, "LL", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);

      }
   }
   
   free(Basis);
   FREE_CRS2(M_CRS, p_threads);
   FREE_CCS2(M_CCS, p_threads);
   FREE_CRS1(M_Base);
   FREE_CRS1(MM_On_T);
   FREE_CRS1(MM_On);
   FREE_ARRAY_DOUBLE2(V, p_threads);
   
}
