#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dmrg.h"

void DMRG_V_M_RL_Q3(CRS1 *M_RL, int qn1_out, int qn2_out, int qn3_out, double *V, int *Ele_LL, int *Ele_LR, int *Ele_RR, char Sign_Flag[], double *Out_V, int p_threads, DMRG_WHOLE_BASIS_Q3 *Dmrg_W_Basis) {
   
   if (qn1_out < 0 || qn2_out < 0 || qn3_out < 0) {
      printf("Error in DMRG_V_M_RL_Q3\n");
      printf("qn1_out=%d,qn2_out=%d,qn3_out=%d\n", qn1_out, qn2_out, qn3_out);
      exit(1);
   }
   
   int dim_RR     = Dmrg_W_Basis->dim_RR;
   int dim_onsite = Dmrg_W_Basis->dim_onsite;
   int basis,LL,LR,RR,RL,col,inv,sign;
   int dim_out = Dmrg_W_Basis->Dim[qn1_out][qn2_out][qn3_out];
   long iter,inv_sup;
   double val;
   
#pragma omp parallel for private (LL,LR,RR,RL,sign,val,iter,col,inv,inv_sup) num_threads (p_threads)
   for (basis = 0; basis < dim_out; basis++) {
      LL = Dmrg_W_Basis->LL_LLLRRRRL[qn1_out][qn2_out][qn3_out][basis];
      LR = Dmrg_W_Basis->LR_LLLRRRRL[qn1_out][qn2_out][qn3_out][basis];
      RR = Dmrg_W_Basis->RR_LLLRRRRL[qn1_out][qn2_out][qn3_out][basis];
      RL = Dmrg_W_Basis->RL_LLLRRRRL[qn1_out][qn2_out][qn3_out][basis];
      
      if (strcmp(Sign_Flag, "Yes") == 0) {
         if ((Ele_LL[LL] + Ele_LR[LR] + Ele_RR[RR] )%2 == 0) {
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
      for (iter = M_RL->Row[RL]; iter < M_RL->Row[RL + 1]; iter++) {
         col     = M_RL->Col[iter];
         inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + col;
         inv     = Dmrg_W_Basis->Inv_LLLRRRRL[inv_sup];
         if (inv >= 0) {
            val = val + sign*V[inv]*M_RL->Val[iter];
         }
      }
      Out_V[basis] = val;
   }
   
}

