#include <stdlib.h>
#include "dmrg.h"

void DMRG_FREE_BASIS_LLLR(DMRG_BASIS *Basis, int dim_LL) {
   
   FREE_ARRAY_SINT1(Basis->LL_LLLR);
   FREE_ARRAY_SINT1(Basis->LR_LLLR);
   FREE_ARRAY_INT2(Basis->Inv_LLLR, dim_LL);
   
   FREE_ARRAY_INT1(Basis->Sum_LLLR);
   
   FREE_ARRAY_INT1(Basis->Tot_Sz_LLLR);
   FREE_ARRAY_INT1(Basis->Tot_Ele_LLLR);
   FREE_ARRAY_INT1(Basis->Tot_Ele_1_LLLR);
   FREE_ARRAY_INT1(Basis->Tot_Ele_2_LLLR);
   FREE_ARRAY_INT1(Basis->Tot_Parity_LLLR);
   
}
