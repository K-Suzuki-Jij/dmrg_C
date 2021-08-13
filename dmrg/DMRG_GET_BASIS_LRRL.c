#include <stdlib.h>
#include <stdio.h>
#include "dmrg.h"

void DMRG_GET_BASIS_LRRL(int dim_onsite, DMRG_BASIS *Basis) {
   
   Basis->LR_LRRL  = GET_ARRAY_SINT1(dim_onsite*dim_onsite);
   Basis->RL_LRRL  = GET_ARRAY_SINT1(dim_onsite*dim_onsite);
   Basis->Inv_LRRL = GET_ARRAY_INT2(dim_onsite, dim_onsite);
   
   int LR,RL;
   for (LR = 0; LR < dim_onsite; LR++) {
      for (RL = 0; RL < dim_onsite; RL++) {
         Basis->Inv_LRRL[LR][RL] = -1;
      }
   }
   
   int basis;
   int dim_LRRL = 0;
   for (basis = 0; basis < Basis->dim_LLLRRRRL; basis++) {
      LR = Basis->LR_LLLRRRRL[basis];
      RL = Basis->RL_LLLRRRRL[basis];
      if (Basis->Inv_LRRL[LR][RL] == -1) {
         Basis->LR_LRRL[dim_LRRL] = LR;
         Basis->RL_LRRL[dim_LRRL] = RL;
         Basis->Inv_LRRL[LR][RL]  = dim_LRRL;
         dim_LRRL++;
      }
   }
   
   if (dim_LRRL <= 0) {
      printf("Error in DMRG_GET_BASIS_LRRL\n");
      printf("dim_LRRL= %d\n",dim_LRRL);
      exit(1);
   }
   
   Basis->dim_LRRL = dim_LRRL;
   
   
}
