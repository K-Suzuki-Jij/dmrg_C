#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "SML.h"
#include "dmrg.h"

void DMRG_TRANS_MAT_C_C(CRS1 *M_On_1, CRS1 *M_On_2, CRS1 **Out, int *Dim_LL, int *Dim_LLLR, CCS1 **TM, CRS1 **TM_D, short **LL_LLLR, short **LR_LLLR, int ***Inv_LLLR, int **Tot_Ele_LL, int tot_site) {
   
   int LL_site,site;
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
   
   CRS1 *M_CRS = GET_CRS1(max_dim_LLLR, max_elem_num_LLLR);
   CCS1 *M_CCS = GET_CCS1(max_dim_LLLR, max_elem_num_LLLR);
   CRS1 *M_LL  = GET_CRS1(max_dim_LL, max_elem_num_LL);
   double *V   = GET_ARRAY_DOUBLE1(max_dim_LLLR);
   
   COPY_CRS1(M_On_1, M_LL, 1);
   DMRG_BASIS_LLLR *Basis = malloc(sizeof(DMRG_BASIS_LLLR));
   
   for (LL_site = 0; LL_site <= tot_site/2 - 3; LL_site++) {
      
      CCS1 *T_M       = TM[LL_site];
      CRS1 *T_MD      = TM_D[LL_site];
      Basis->dim_LLLR = Dim_LLLR[LL_site];
      Basis->LL_LLLR  = LL_LLLR[LL_site];
      Basis->LR_LLLR  = LR_LLLR[LL_site];
      Basis->Inv_LLLR = Inv_LLLR[LL_site];
      
      DMRG_EXTEND_AND_REDUCE(M_LL, M_On_2, Out[LL_site], Tot_Ele_LL[LL_site], "LL_LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Basis);
      
      for (site = 0; site < LL_site; site++) {
         DMRG_EXTEND_AND_REDUCE(Out[site], NULL, Out[site], NULL, "LL", "No", T_M, T_MD, M_CRS, M_CCS, V, Basis);
      }
      
      DMRG_EXTEND_AND_REDUCE(NULL, M_On_1, M_LL, Tot_Ele_LL[LL_site], "LR", "Yes", T_M, T_MD, M_CRS, M_CCS, V, Basis);
      
   }
   
   FREE_CRS1(M_CRS);
   FREE_CCS1(M_CCS);
   FREE_CRS1(M_LL);
   FREE_ARRAY_DOUBLE1(V);
   free(Basis);
   
}
