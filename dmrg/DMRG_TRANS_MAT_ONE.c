#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SML.h"
#include "dmrg.h"

void DMRG_TRANS_MAT_ONE(CRS1 *M_On, CRS1 **Out, int *Dim_LL, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int tot_site, int p_threads) {
   
   int site,thread_num,LL_site;
   int max_elem_num_LLLR = 0;
   int max_dim_LLLR      = 0;

   for (LL_site = 0; LL_site <= tot_site/2 - 2; LL_site++) {
      if (max_elem_num_LLLR < Dim_LLLR[LL_site]*Dim_LL[LL_site + 1]) {
         max_elem_num_LLLR = Dim_LLLR[LL_site]*Dim_LL[LL_site + 1];
      }
      if (max_dim_LLLR < Dim_LLLR[LL_site]) {
         max_dim_LLLR = Dim_LLLR[LL_site];
      }
   }
   
   CRS1 **M_CRS = GET_CRS2(p_threads, max_dim_LLLR, max_elem_num_LLLR);
   CCS1 **M_CCS = GET_CCS2(p_threads, max_dim_LLLR, max_elem_num_LLLR);
   double **V   = GET_ARRAY_DOUBLE2(p_threads, max_dim_LLLR);
   
   COPY_CRS1(M_On, Out[0], 1);
   
   DMRG_BASIS_LLLR *Basis = malloc(sizeof(DMRG_BASIS_LLLR));
   
   for (LL_site = 0; LL_site <= tot_site/2 - 3; LL_site++) {
      
      CCS1 *T_M       = TM[LL_site];
      CRS1 *T_MD      = TM_D[LL_site];
      Basis->dim_LLLR = Dim_LLLR[LL_site];
      Basis->LL_LLLR  = LL_LLLR[LL_site];
      Basis->LR_LLLR  = LR_LLLR[LL_site];
      Basis->Inv_LLLR = Inv_LLLR[LL_site];
      
#pragma omp parallel for private (thread_num) num_threads (p_threads)
      for (site = 0; site <= LL_site; site++) {
         thread_num = omp_get_thread_num();
         DMRG_EXTEND_AND_REDUCE(Out[site], NULL, Out[site], NULL, "LL", "No", T_M, T_MD, M_CRS[thread_num], M_CCS[thread_num], V[thread_num], Basis);
      }
      
      DMRG_EXTEND_AND_REDUCE(NULL, M_On, Out[LL_site + 1], NULL, "LR", "No", T_M, T_MD, M_CRS[0], M_CCS[0], V[0], Basis);
      
   }
   
   free(Basis);
   FREE_CRS2(M_CRS, p_threads);
   FREE_CCS2(M_CCS, p_threads);
   FREE_ARRAY_DOUBLE2(V, p_threads);
   
}
