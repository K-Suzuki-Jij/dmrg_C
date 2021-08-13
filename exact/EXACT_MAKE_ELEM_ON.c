#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SML.h"
#include "exact.h"

void EXACT_MAKE_ELEM_ON(long basis, int site, int dim_onsite, CRS1 *M_On, long *elem_num, double coeef, EXACT_A_BASIS *A_Basis) {

   int local_basis = EXACT_FIND_SITE_STATE(basis, site, dim_onsite);
   long dim_site   = (long)pow(dim_onsite, site);
   long temp_elem_num = *elem_num;
   long whole_a_basis,iter,inv;
   
   for (iter = M_On->Row[local_basis]; iter < M_On->Row[local_basis + 1]; iter++) {
      whole_a_basis = basis + (M_On->Col[iter] - local_basis)*dim_site;
      inv           = BINARY_SEARCH_LINT1(A_Basis->Check, 0, temp_elem_num, whole_a_basis);
      if (inv == -1) {
         A_Basis->Basis[temp_elem_num] = whole_a_basis;
         A_Basis->Check[temp_elem_num] = whole_a_basis;
         A_Basis->Val[temp_elem_num]   = coeef*M_On->Val[iter];
         temp_elem_num++;
         QUICK_SORT_LINT2_DOUBLE1(A_Basis->Check, A_Basis->Basis, A_Basis->Val, 0, temp_elem_num);
      }
      else {
         A_Basis->Val[inv] = A_Basis->Val[inv] + coeef*M_On->Val[iter];
      }
   }
   
   *elem_num = temp_elem_num;
   
}
