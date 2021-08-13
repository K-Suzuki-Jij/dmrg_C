//
//  GET_HAM_LLLRRRRL.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

CRS1 *GET_HAM_LLLRRRRL(DMRG_BASIS *Dmrg_Basis, BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->make_ham = omp_get_wtime();
   
   int dim_RR       = Dmrg_Basis->dim_RR;
   int dim_onsite   = Dmrg_Basis->dim_onsite;
   int LL_site      = Dmrg_Status->LL_site;
   int RR_site      = Dmrg_Status->RR_site;
   int dim_LLLRRRRL = Dmrg_Basis->dim_LLLRRRRL;

   DMRG_A_BASIS **Dmrg_A_Basis = DMRG_GET_A_BASIS(dim_LLLRRRRL, Model->p_threads);
   HAM_BOX  **Box              = GET_HAM_BOX(System, Enviro, Model, LL_site, RR_site);
   
   long *Row_Elem_Num = GET_ARRAY_LINT1(dim_LLLRRRRL + 1);
   int basis,thread_num,num,inv,LL,LR,RR,RL;
   double zero = pow(10,-15);
   double val;
   long inv_sup;
   DMRG_BASIS_ONSITE Basis;

#pragma omp parallel for private (thread_num,Basis,num,LL,LR,RR,RL,inv,val,inv_sup) schedule(guided) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      thread_num       = omp_get_thread_num();
      Basis.LL         = Dmrg_Basis->LL_LLLRRRRL[basis];
      Basis.LR         = Dmrg_Basis->LR_LLLRRRRL[basis];
      Basis.RL         = Dmrg_Basis->RL_LLLRRRRL[basis];
      Basis.RR         = Dmrg_Basis->RR_LLLRRRRL[basis];
      Basis.dim_RR     = dim_RR;
      Basis.dim_onsite = dim_onsite;
      MAKE_ELEMENT_HAM_LLLRRRRL(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_LLLRRRRL, Box[thread_num]);
      for (num = 0; num < Dmrg_A_Basis[thread_num]->elem_num; num++) {
         LL      = Dmrg_A_Basis[thread_num]->LL_LLLRRRRL[num];
         LR      = Dmrg_A_Basis[thread_num]->LR_LLLRRRRL[num];
         RR      = Dmrg_A_Basis[thread_num]->RR_LLLRRRRL[num];
         RL      = Dmrg_A_Basis[thread_num]->RL_LLLRRRRL[num];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         val     = Dmrg_A_Basis[thread_num]->Val_LLLRRRRL[num];
         if ((0 <= inv && inv <= basis && fabs(val) > zero) || (inv == basis)) {
            Row_Elem_Num[basis + 1]++;
         }
      }
   }
  
   long tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLRRRRL + 1; basis++) {
      tot_elem_num = tot_elem_num + Row_Elem_Num[basis];
   }
   
   CRS1 *Ham = GET_CRS1(dim_LLLRRRRL, tot_elem_num + 1);

   //Do not use openmp here
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      Row_Elem_Num[basis + 1] = Row_Elem_Num[basis + 1] + Row_Elem_Num[basis];
   }
   
#pragma omp parallel for private (thread_num,Basis,num,LL,LR,RR,RL,inv,val,inv_sup) schedule(guided) num_threads (Model->p_threads)
   for (basis = 0; basis < dim_LLLRRRRL; basis++) {
      thread_num       = omp_get_thread_num();
      Basis.LL         = Dmrg_Basis->LL_LLLRRRRL[basis];
      Basis.LR         = Dmrg_Basis->LR_LLLRRRRL[basis];
      Basis.RL         = Dmrg_Basis->RL_LLLRRRRL[basis];
      Basis.RR         = Dmrg_Basis->RR_LLLRRRRL[basis];
      Basis.dim_RR     = dim_RR;
      Basis.dim_onsite = dim_onsite;
      MAKE_ELEMENT_HAM_LLLRRRRL(&Basis, Dmrg_A_Basis[thread_num], Dmrg_Basis->Inv_LLLRRRRL, Box[thread_num]);
      for (num = 0; num < Dmrg_A_Basis[thread_num]->elem_num; num++) {
         LL      = Dmrg_A_Basis[thread_num]->LL_LLLRRRRL[num];
         LR      = Dmrg_A_Basis[thread_num]->LR_LLLRRRRL[num];
         RR      = Dmrg_A_Basis[thread_num]->RR_LLLRRRRL[num];
         RL      = Dmrg_A_Basis[thread_num]->RL_LLLRRRRL[num];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
         inv     = Dmrg_Basis->Inv_LLLRRRRL[inv_sup];
         val     = Dmrg_A_Basis[thread_num]->Val_LLLRRRRL[num];
         if ((0 <= inv && inv <= basis && fabs(val) > zero) || (inv == basis)) {
            Ham->Col[Row_Elem_Num[basis]] = inv;
            Ham->Val[Row_Elem_Num[basis]] = val;
            Row_Elem_Num[basis] = Row_Elem_Num[basis] + 1;
         }
      }
      Ham->Row[basis + 1] = Row_Elem_Num[basis];
   }
   
   if (Ham->Row[dim_LLLRRRRL] != tot_elem_num) {
      printf("Error in GET_HAM_LLLRRRRL\n");
      printf("%ld != %ld\n", Ham->Row[dim_LLLRRRRL], tot_elem_num);
      exit(1);
   }

   Ham->row_dim = dim_LLLRRRRL;
   Ham->col_dim = dim_LLLRRRRL;
   
   SORT_COLUMN_CRS1(Ham, Model->p_threads);
   
   DMRG_FREE_A_BASIS(Dmrg_A_Basis, Model->p_threads);
   FREE_ARRAY_LINT1(Row_Elem_Num);
   FREE_HAM_BOX(Box, Model);
   
   Dmrg_Time->make_ham = omp_get_wtime() - Dmrg_Time->make_ham;
   return Ham;
}

