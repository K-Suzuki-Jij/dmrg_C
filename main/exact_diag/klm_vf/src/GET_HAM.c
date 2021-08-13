//
//  GET_HAM.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

CRS1 *GET_HAM(MODEL_1DKLM_VF *Model, EXACT_BASIS_INFO *Basis_Info, EXACT_PARAMETER *Param, EXACT_TIME *Time) {
   
   Time->make_ham = omp_get_wtime();
   
   HAM_BOX *Ham_Box        = GET_HAM_BOX(Model);
   EXACT_A_BASIS **A_Basis = EXACT_GET_A_BASIS(Model->p_threads, Param->est_max_row_elem_num);
   int dim                 = Basis_Info->dim;
   long *Row_Elem_Num      = GET_ARRAY_LINT1(dim + 1);
   double zero             = pow(10,-15);
   int basis,thread_num,i;
   long tot_elem_num,inv;
   double val;
   
#pragma omp parallel for private (thread_num,i,inv,val) num_threads (Model->p_threads)
   for (basis = 0; basis < dim; basis++) {
      thread_num = omp_get_thread_num();
      MAKE_ELEMENT_HAM(Basis_Info->Basis[basis], A_Basis[thread_num], Ham_Box, Model);
      if (Param->est_max_row_elem_num <= A_Basis[thread_num]->elem_num) {
         printf("Error in GET_HAM\n");
         exit(1);
      }
      for (i = 0; i < A_Basis[thread_num]->elem_num; i++) {
         inv = BINARY_SEARCH_LINT1(Basis_Info->Basis, 0, dim, A_Basis[thread_num]->Basis[i]);
         val = A_Basis[thread_num]->Val[i];
         if ((fabs(val) > zero && inv >= 0) || inv == basis) {
            Row_Elem_Num[basis + 1]++;
         }
      }
   }
   
   tot_elem_num = 0;
   
#pragma omp parallel for reduction(+:tot_elem_num) num_threads (Model->p_threads)
   for (basis = 0; basis < dim + 1; basis++) {
      tot_elem_num = tot_elem_num + Row_Elem_Num[basis];
   }
   
   //Do not use openmp here
   for (basis = 0; basis < dim; basis++) {
      Row_Elem_Num[basis + 1] = Row_Elem_Num[basis + 1] + Row_Elem_Num[basis];
   }
   
   //Make Hamiltonian
   CRS1 *Out = GET_CRS1(dim, tot_elem_num);
   
#pragma omp parallel for private (thread_num,i,inv,val) num_threads (Model->p_threads)
   for (basis = 0; basis < dim; basis++) {
      thread_num = omp_get_thread_num();
      MAKE_ELEMENT_HAM(Basis_Info->Basis[basis], A_Basis[thread_num], Ham_Box, Model);
      for (i = 0; i < A_Basis[thread_num]->elem_num; i++) {
         inv = BINARY_SEARCH_LINT1(Basis_Info->Basis, 0, dim, A_Basis[thread_num]->Basis[i]);
         val = A_Basis[thread_num]->Val[i];
         if ((fabs(val) > zero && inv >= 0) || inv == basis) {
            Out->Val[Row_Elem_Num[basis]] = val;
            Out->Col[Row_Elem_Num[basis]] = (int)inv;
            Row_Elem_Num[basis] = Row_Elem_Num[basis] + 1;
         }
      }
      Out->Row[basis + 1] = Row_Elem_Num[basis];
   }
   
   if (Out->Row[dim] != tot_elem_num) {
      printf("Error in GET_HAM\n");
      printf("%ld != %ld\n",Out->Row[dim], tot_elem_num);
      exit(1);
   }
   
   Out->row_dim = dim;
   Out->col_dim = dim;
   
   SORT_COLUMN_CRS1(Out, Model->p_threads);
   FREE_ARRAY_LINT1(Row_Elem_Num);
   FREE_HAM_BOX(Ham_Box, Model);
   EXACT_FREE_A_BASIS(A_Basis, Model->p_threads);
   
   CHECK_SYMMETRY_CRS1(Out, pow(10,-15), 1);
   
   Time->make_ham = omp_get_wtime() - Time->make_ham;
   
   return Out;
   
}
