#include <stdlib.h>
#include "dmrg.h"

void DMRG_GET_BASIS_LLLRRRRL(long dim_LLLRRRRL, int dim_onsite, int dim_LL, int dim_RR, DMRG_BASIS *Basis) {
   
   Basis->LL_LLLRRRRL  = GET_ARRAY_SINT1(dim_LLLRRRRL + 10);
   Basis->LR_LLLRRRRL  = GET_ARRAY_SINT1(dim_LLLRRRRL + 10);
   Basis->RL_LLLRRRRL  = GET_ARRAY_SINT1(dim_LLLRRRRL + 10);
   Basis->RR_LLLRRRRL  = GET_ARRAY_SINT1(dim_LLLRRRRL + 10);
   Basis->Inv_LLLRRRRL = GET_ARRAY_INT1((long)dim_LL*dim_onsite*dim_RR*dim_onsite);
   
}
