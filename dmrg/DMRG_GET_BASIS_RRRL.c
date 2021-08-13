#include <stdlib.h>
#include <stdio.h>
#include "dmrg.h"

void DMRG_GET_BASIS_RRRL(int dim_RR, int dim_onsite, DMRG_BASIS *Basis) {
   
   Basis->RR_RRRL  = GET_ARRAY_SINT1(dim_RR*dim_onsite);
   Basis->RL_RRRL  = GET_ARRAY_SINT1(dim_RR*dim_onsite);
   Basis->Inv_RRRL = GET_ARRAY_INT2(dim_RR, dim_onsite);
   
   int RR,RL;
   for (RR = 0; RR < dim_RR; RR++) {
      for (RL = 0; RL < dim_onsite; RL++) {
         Basis->Inv_RRRL[RR][RL] = -1;
      }
   }
   
   int basis;
   int dim_RRRL = 0;
   for (basis = 0; basis < Basis->dim_LLLRRRRL; basis++) {
      RR = Basis->RR_LLLRRRRL[basis];
      RL = Basis->RL_LLLRRRRL[basis];
      if (Basis->Inv_RRRL[RR][RL] == -1) {
         Basis->RR_RRRL[dim_RRRL] = RR;
         Basis->RL_RRRL[dim_RRRL] = RL;
         Basis->Inv_RRRL[RR][RL]  = dim_RRRL;
         dim_RRRL++;
      }
   }
   
   if (dim_RRRL <= 0) {
      printf("Error in DMRG_GET_BASIS_RRRL\n");
      printf("dim_RRRL= %d\n",dim_RRRL);
      exit(1);
   }
   
   Basis->dim_RRRL = dim_RRRL;
   
}
