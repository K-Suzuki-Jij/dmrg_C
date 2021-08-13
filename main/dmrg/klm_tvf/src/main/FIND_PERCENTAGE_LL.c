//
//  FIND_PERCENTAGE_LL.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DKLM_TVF *Model) {
   
   int tot_site  = Model->tot_site;
   int tot_num   = 30*tot_site;
   double *P     = GET_ARRAY_DOUBLE1(tot_num);
   int site;
   
   for (site = 0; site < tot_site; site++) {
      P[site + 0*tot_site] = (double)System->Ham[site]->Row[System->Ham[site]->row_dim]/System->Ham[site]->max_val;
      P[site + 1*tot_site] = (double)System->SzL_RE[site]->Row[System->SzL_RE[site]->row_dim]/System->SzL_RE[site]->max_val;
      P[site + 2*tot_site] = (double)System->SpL_RE[site]->Row[System->SpL_RE[site]->row_dim]/System->SpL_RE[site]->max_val;
      P[site + 4*tot_site] = (double)System->SmL_RE[site]->Row[System->SmL_RE[site]->row_dim]/System->SmL_RE[site]->max_val;
      P[site + 5*tot_site] = (double)System->Even_RE[site]->Row[System->Even_RE[site]->row_dim]/System->Even_RE[site]->max_val;
      P[site + 6*tot_site] = (double)System->Odd_RE[site]->Row[System->Odd_RE[site]->row_dim]/System->Odd_RE[site]->max_val;
      P[site + 7*tot_site] = (double)System->Even_D_RE[site]->Row[System->Even_D_RE[site]->row_dim]/System->Even_D_RE[site]->max_val;
      P[site + 8*tot_site] = (double)System->Odd_D_RE[site]->Row[System->Odd_D_RE[site]->row_dim]/System->Odd_D_RE[site]->max_val;
   }
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      for (site = 0; site < tot_site; site++) {
         P[site + 9*tot_site]  = (double)System->SzL_LE[site]->Row[System->SzL_LE[site]->row_dim]/System->SzL_LE[site]->max_val;
         P[site + 10*tot_site] = (double)System->SpL_LE[site]->Row[System->SpL_LE[site]->row_dim]/System->SpL_LE[site]->max_val;
         P[site + 11*tot_site] = (double)System->SmL_LE[site]->Row[System->SmL_LE[site]->row_dim]/System->SmL_LE[site]->max_val;
         P[site + 12*tot_site] = (double)System->Even_LE[site]->Row[System->Even_LE[site]->row_dim]/System->Even_LE[site]->max_val;
         P[site + 13*tot_site] = (double)System->Odd_LE[site]->Row[System->Odd_LE[site]->row_dim]/System->Odd_LE[site]->max_val;
         P[site + 14*tot_site] = (double)System->Even_D_LE[site]->Row[System->Even_D_LE[site]->row_dim]/System->Even_D_LE[site]->max_val;
         P[site + 15*tot_site] = (double)System->Odd_D_LE[site]->Row[System->Odd_D_LE[site]->row_dim]/System->Odd_D_LE[site]->max_val;
      }
   }
   
   double max = 0;
   for (site = 0; site < tot_num; site++) {
      if (max < P[site]) {
         max = P[site];
      }
   }

   FREE_ARRAY_DOUBLE1(P);
   
   return max;
   
}
