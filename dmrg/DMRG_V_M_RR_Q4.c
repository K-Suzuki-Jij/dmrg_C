#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dmrg.h"

void DMRG_V_M_RR_Q4(CRS1 *M_RR, int qn1_out, int qn2_out, int qn3_out, int qn4_out, double *V, int *Ele_LL, int *Ele_LR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q4 *Dmrg_W_Basis) {
   
   if (qn1_out < 0 || qn2_out < 0 || qn3_out < 0 || qn4_out < 0) {
      printf("Error in DMRG_V_M_RR_Q4\n");
      printf("qn1_out=%d,qn2_out=%d,qn3_out=%d,qn4_out=%d\n", qn1_out, qn2_out, qn3_out, qn4_out);
      exit(1);
   }
   int dim_RR     = Dmrg_W_Basis->dim_RR;
   int dim_onsite = Dmrg_W_Basis->dim_onsite;
   int basis,LL,LR,RR,RL,col,inv,sign;
   int dim_out = Dmrg_W_Basis->Dim[qn1_out][qn2_out][qn3_out][qn4_out];
   long iter,inv_sup;
   double val;
   
#pragma omp parallel for private (LL,LR,RR,RL,sign,val,iter,col,inv,inv_sup) num_threads (p_threads)
   for (basis = 0; basis < dim_out; basis++) {
      LL = Dmrg_W_Basis->LL_LLLRRRRL[qn1_out][qn2_out][qn3_out][qn4_out][basis];
      LR = Dmrg_W_Basis->LR_LLLRRRRL[qn1_out][qn2_out][qn3_out][qn4_out][basis];
      RR = Dmrg_W_Basis->RR_LLLRRRRL[qn1_out][qn2_out][qn3_out][qn4_out][basis];
      RL = Dmrg_W_Basis->RL_LLLRRRRL[qn1_out][qn2_out][qn3_out][qn4_out][basis];
      
      if (strcmp(Sign_Flag, "Yes") == 0) {
         if ((Ele_LL[LL] + Ele_LR[LR])%2 == 0) {
            sign = 1;
         }
         else {
            sign = -1;
         }
      }
      else {
         sign = 1;
      }
      
      val = 0;
      for (iter = M_RR->Row[RR]; iter < M_RR->Row[RR + 1]; iter++) {
         col     = M_RR->Col[iter];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + col*dim_onsite + RL;
         inv     = Dmrg_W_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + sign*V[inv]*M_RR->Val[iter];
         }
      }
      Out_V[basis] = val;
   }
   
}

