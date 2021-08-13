//
//  GET_HAM_LLLR.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

CRS1 *GET_HAM_LLLR(BLOCK *System, BLOCK *Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DHUBBARD_VF *Model) {
   
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int dim_LLLR   = Dmrg_Basis->dim_LLLR;
   
   
   DMRG_A_BASIS **Dmrg_A_Basis = DMRG_GET_A_BASIS(dim_LLLR, Model->p_threads);
   HAM_BOX  **Box              = GET_HAM_BOX(System, Enviro, Model, LL_site, RR_site);
   
   long *Row_Elem_Num = GET_ARRAY_LINT1(dim_LLLR + 1);
   int num,inv,LL,LR,basis,thread_num;
   double zero = pow(10,-15);
   double val;
   DMRG_BASIS_ONSITE Basis;
   
#pragma omp parallel for private (thread_num,Basis,num,LL,LR,inv,val) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLR; basis++) {
      thread_num = omp_get_thread_num();
      Basis.LL   = Dmrg_Basis->LL_LLLR[basis];
      Basis.LR   = Dmrg_Basis->LR_LLLR[basis];
      MAKE_ELEMENT_HAM_LLLR(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_LLLR, Box[thread_num]);
      for (num = 0; num < Dmrg_A_Basis[thread_num]->elem_num; num++) {
         LL  = Dmrg_A_Basis[thread_num]->LL_LLLRRRRL[num];
         LR  = Dmrg_A_Basis[thread_num]->LR_LLLRRRRL[num];
         inv = Dmrg_Basis->Inv_LLLR[LL][LR];
         val = Dmrg_A_Basis[thread_num]->Val_LLLRRRRL[num];
         if (inv >= 0 && fabs(val) > zero) {
            Row_Elem_Num[basis + 1]++;
         }
      }
   }
   
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLR + 1; basis++) {
      tot_elem_num = tot_elem_num + Row_Elem_Num[basis];
   }
   
   //Do not use openmp here
   for (basis = 0; basis < dim_LLLR; basis++) {
      Row_Elem_Num[basis + 1] = Row_Elem_Num[basis + 1] + Row_Elem_Num[basis];
   }
   
   CRS1 *Ham_LLLR = GET_CRS1(dim_LLLR, tot_elem_num);
   
#pragma omp parallel for private (thread_num,Basis,num,LL,LR,inv,val) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLR; basis++) {
      thread_num = omp_get_thread_num();
      Basis.LL   = Dmrg_Basis->LL_LLLR[basis];
      Basis.LR   = Dmrg_Basis->LR_LLLR[basis];
      MAKE_ELEMENT_HAM_LLLR(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_LLLR, Box[thread_num]);
      for (num = 0; num < Dmrg_A_Basis[thread_num]->elem_num; num++) {
         LL  = Dmrg_A_Basis[thread_num]->LL_LLLRRRRL[num];
         LR  = Dmrg_A_Basis[thread_num]->LR_LLLRRRRL[num];
         inv = Dmrg_Basis->Inv_LLLR[LL][LR];
         val = Dmrg_A_Basis[thread_num]->Val_LLLRRRRL[num];
         if (inv >= 0 && fabs(val) > zero) {
            Ham_LLLR->Col[Row_Elem_Num[basis]] = inv;
            Ham_LLLR->Val[Row_Elem_Num[basis]] = val;
            Row_Elem_Num[basis] = Row_Elem_Num[basis] + 1;
         }
      }
      Ham_LLLR->Row[basis + 1] = Row_Elem_Num[basis];
   }
   
   //Check Point
   if (Ham_LLLR->Row[dim_LLLR] != tot_elem_num) {
      printf("Error in GET_HAM_LLLR\n");
      printf("%ld != %ld\n", Ham_LLLR->Row[dim_LLLR], tot_elem_num);
      exit(1);
   }
   
   Ham_LLLR->row_dim = dim_LLLR;
   Ham_LLLR->col_dim = dim_LLLR;
   
   SORT_COLUMN_CRS1(Ham_LLLR, Model->p_threads);
   
   DMRG_FREE_A_BASIS(Dmrg_A_Basis, Model->p_threads);
   FREE_HAM_BOX(Box, Model);
   FREE_ARRAY_LINT1(Row_Elem_Num);
   
   return Ham_LLLR;
}
