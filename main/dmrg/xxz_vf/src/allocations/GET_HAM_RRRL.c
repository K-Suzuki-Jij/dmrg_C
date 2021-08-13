//
//  GET_HAM_RRRL.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/27.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

CRS1 *GET_HAM_RRRL(BLOCK *Block_System, BLOCK *Block_Enviro, DMRG_BASIS *Dmrg_Basis, DMRG_STATUS *Dmrg_Status, MODEL_1DXXZ_VF *Model) {
   
   int LL_site    = Dmrg_Status->LL_site;
   int RR_site    = Dmrg_Status->RR_site;
   int dim_RRRL   = Dmrg_Basis->dim_RRRL;
   
   DMRG_A_BASIS **Dmrg_A_Basis = DMRG_GET_A_BASIS(dim_RRRL, Model->p_threads);
   HAM_BOX  **Box              = GET_HAM_BOX(Block_System, Block_Enviro, Model, LL_site, RR_site);
   
   long *Row_Elem_Num = GET_ARRAY_LINT1(dim_RRRL + 1);
   int num,inv,RR,RL,basis,thread_num;
   DMRG_BASIS_ONSITE Basis;
   
#pragma omp parallel for private (thread_num,Basis) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_RRRL; basis++) {
      thread_num = omp_get_thread_num();
      Basis.RR   = Dmrg_Basis->RR_RRRL[basis];
      Basis.RL   = Dmrg_Basis->RL_RRRL[basis];
      MAKE_ELEMENT_HAM_RRRL(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_RRRL, Box[thread_num]);
      Row_Elem_Num[basis + 1] = Dmrg_A_Basis[thread_num]->elem_num;
   }

   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_RRRL + 1; basis++) {
      tot_elem_num = tot_elem_num + Row_Elem_Num[basis];
   }
   
   //Do not use openmp here
   for (basis = 0; basis < dim_RRRL; basis++) {
      Row_Elem_Num[basis + 1] = Row_Elem_Num[basis + 1] + Row_Elem_Num[basis];
   }
   
   CRS1 *Ham_RRRL = GET_CRS1(dim_RRRL, tot_elem_num);
   
#pragma omp parallel for private (thread_num,Basis,num,RR,RL,inv) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_RRRL; basis++) {
      thread_num = omp_get_thread_num();
      Basis.RR   = Dmrg_Basis->RR_RRRL[basis];
      Basis.RL   = Dmrg_Basis->RL_RRRL[basis];
      MAKE_ELEMENT_HAM_RRRL(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_RRRL, Box[thread_num]);
      for (num = 0; num < Dmrg_A_Basis[thread_num]->elem_num; num++) {
         RR  = Dmrg_A_Basis[thread_num]->RR_LLLRRRRL[num];
         RL  = Dmrg_A_Basis[thread_num]->RL_LLLRRRRL[num];
         inv = Dmrg_Basis->Inv_RRRL[RR][RL];
         if (inv >= 0) {
            Ham_RRRL->Col[Row_Elem_Num[basis]] = inv;
            Ham_RRRL->Val[Row_Elem_Num[basis]] = Dmrg_A_Basis[thread_num]->Val_LLLRRRRL[num];
            Row_Elem_Num[basis] = Row_Elem_Num[basis] + 1;
         }
         else {
            printf("Error in GET_HAM_RRRL\n");
            printf("inv = %d\n",inv);
            exit(1);
         }
      }
      Ham_RRRL->Row[basis + 1] = Dmrg_A_Basis[thread_num]->elem_num;
   }
      
   //Do not use openmp here
   for (basis = 0; basis < dim_RRRL; basis++) {
      Ham_RRRL->Row[basis + 1] = Ham_RRRL->Row[basis + 1] + Ham_RRRL->Row[basis];
   }
   
   //Check Point
   if (Ham_RRRL->Row[dim_RRRL] != tot_elem_num) {
      printf("Error in GET_HAM_RRRL\n");
      printf("%ld != %ld\n", Ham_RRRL->Row[dim_RRRL], tot_elem_num);
      exit(1);
   }
   
   Ham_RRRL->row_dim = dim_RRRL;
   Ham_RRRL->col_dim = dim_RRRL;
   
   SORT_COLUMN_CRS1(Ham_RRRL, Model->p_threads);
   
   DMRG_FREE_A_BASIS(Dmrg_A_Basis, Model->p_threads);
   FREE_HAM_BOX(Box, Model);
   FREE_ARRAY_LINT1(Row_Elem_Num);
   
   return Ham_RRRL;
}
