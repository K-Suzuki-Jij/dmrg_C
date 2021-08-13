#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "exact.h"

void EXACT_MAKE_ELEM_INTER(long basis, int site1, int site2, int dim_onsite, CRS1 *M1, CRS1 *M2, long *elem_num, double coeef, int sign, EXACT_A_BASIS *A_Basis) {
   
   int local_basis1   = EXACT_FIND_SITE_STATE(basis, site1, dim_onsite);
   int local_basis2   = EXACT_FIND_SITE_STATE(basis, site2, dim_onsite);
   long dim_site1     = (long)pow(dim_onsite, site1);
   long dim_site2     = (long)pow(dim_onsite, site2);
   long temp_elem_num = *elem_num;
   long whole_a_basis,iter1,iter2,inv;
   double val;
   
   for (iter1 = M1->Row[local_basis1]; iter1 < M1->Row[local_basis1 + 1]; iter1++) {
      for (iter2 = M2->Row[local_basis2]; iter2 < M2->Row[local_basis2 + 1]; iter2++) {
         whole_a_basis   = basis + (M1->Col[iter1] - local_basis1)*dim_site1 + (M2->Col[iter2] - local_basis2)*dim_site2;
         inv             = BINARY_SEARCH_LINT1(A_Basis->Check, 0, temp_elem_num, whole_a_basis);
         val             = sign*coeef*M1->Val[iter1]*M2->Val[iter2];
         if (inv == -1) {
            A_Basis->Basis[temp_elem_num] = whole_a_basis;
            A_Basis->Check[temp_elem_num] = whole_a_basis;
            A_Basis->Val[temp_elem_num]   = val;
            temp_elem_num++;
            QUICK_SORT_LINT2_DOUBLE1(A_Basis->Check, A_Basis->Basis, A_Basis->Val, 0, temp_elem_num);
         }
         else {
            A_Basis->Val[inv] = A_Basis->Val[inv] + val;
         }
      }
   }
   
   *elem_num = temp_elem_num;
   
}
