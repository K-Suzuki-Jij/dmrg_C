#include <stdlib.h>
#include "dmrg.h"

void DMRG_GET_BASIS_LLLR(int dim_onsite, int dim_LL, DMRG_BASIS *Basis) {
   
   Basis->LL_LLLR  = GET_ARRAY_SINT1(dim_LL*dim_onsite);
   Basis->LR_LLLR  = GET_ARRAY_SINT1(dim_LL*dim_onsite);
   Basis->Inv_LLLR = GET_ARRAY_INT2(dim_LL, dim_onsite);
   
   Basis->Sum_LLLR = GET_ARRAY_INT1(dim_LL*dim_onsite);
   
   Basis->Tot_Sz_LLLR     = GET_ARRAY_INT1(dim_LL*dim_onsite);
   Basis->Tot_Ele_LLLR    = GET_ARRAY_INT1(dim_LL*dim_onsite);
   Basis->Tot_Ele_1_LLLR  = GET_ARRAY_INT1(dim_LL*dim_onsite);
   Basis->Tot_Ele_2_LLLR  = GET_ARRAY_INT1(dim_LL*dim_onsite);
   Basis->Tot_Parity_LLLR = GET_ARRAY_INT1(dim_LL*dim_onsite);
   
   int LL,LR;
   for (LL = 0; LL < dim_LL; LL++) {
      for (LR = 0; LR < dim_onsite; LR++) {
         Basis->Inv_LLLR[LL][LR] = -1;
      }
   }
   
}
