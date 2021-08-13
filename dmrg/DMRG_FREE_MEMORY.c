#include "dmrg.h"
#include <stdlib.h>

void DMRG_FREE_MEMORY(DMRG_BASIS *Dmrg_Basis, DMRG_SYSTEM_INFO *Dmrg_System) {
   
   DMRG_FREE_BASIS_LLLRRRRL(Dmrg_Basis);
   DMRG_FREE_BASIS_LLLR(Dmrg_Basis, Dmrg_Basis->dim_LL);
   FREE_ARRAY_SINT1(Dmrg_Basis->LR_LRRL);
   FREE_ARRAY_SINT1(Dmrg_Basis->RL_LRRL);
   FREE_ARRAY_SINT1(Dmrg_Basis->RR_RRRL);
   FREE_ARRAY_SINT1(Dmrg_Basis->RL_RRRL);
   FREE_ARRAY_INT2(Dmrg_Basis->Inv_LRRL, Dmrg_Basis->dim_onsite);
   FREE_ARRAY_INT2(Dmrg_Basis->Inv_RRRL, Dmrg_Basis->dim_RR);
   free(Dmrg_Basis);
   DMRG_FREE_SYSTEM_INFO(Dmrg_System);
   
   
}
