//
//  FIND_PERCENTAGE_LL.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/22.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DTKLM_VF *Model) {
   
   int tot_site  = Model->tot_site;
   int tot_num   = 30*tot_site;
   double *P     = GET_ARRAY_DOUBLE1(tot_num);
   int site;
   
   for (site = 0; site < tot_site; site++) {
      P[site + 0*tot_site ] = (double)System->Ham[site]->Row[System->Ham[site]->row_dim]/System->Ham[site]->max_val;
      P[site + 1*tot_site ] = (double)System->SzL_RE[site]->Row[System->SzL_RE[site]->row_dim]/System->SzL_RE[site]->max_val;
      P[site + 2*tot_site ] = (double)System->SpL_RE[site]->Row[System->SpL_RE[site]->row_dim]/System->SpL_RE[site]->max_val;
      P[site + 4*tot_site ] = (double)System->SmL_RE[site]->Row[System->SmL_RE[site]->row_dim]/System->SmL_RE[site]->max_val;
      P[site + 5*tot_site ] = (double)System->CUp_1_RE[site]->Row[System->CUp_1_RE[site]->row_dim]/System->CUp_1_RE[site]->max_val;
      P[site + 6*tot_site ] = (double)System->CDown_1_RE[site]->Row[System->CDown_1_RE[site]->row_dim]/System->CDown_1_RE[site]->max_val;
      P[site + 7*tot_site ] = (double)System->CUp_1_D_RE[site]->Row[System->CUp_1_D_RE[site]->row_dim]/System->CUp_1_D_RE[site]->max_val;
      P[site + 8*tot_site ] = (double)System->CDown_1_D_RE[site]->Row[System->CDown_1_D_RE[site]->row_dim]/System->CDown_1_D_RE[site]->max_val;
      P[site + 9*tot_site ] = (double)System->CUp_2_RE[site]->Row[System->CUp_2_RE[site]->row_dim]/System->CUp_2_RE[site]->max_val;
      P[site + 10*tot_site] = (double)System->CDown_2_RE[site]->Row[System->CDown_2_RE[site]->row_dim]/System->CDown_2_RE[site]->max_val;
      P[site + 11*tot_site] = (double)System->CUp_2_D_RE[site]->Row[System->CUp_2_D_RE[site]->row_dim]/System->CUp_2_D_RE[site]->max_val;
      P[site + 12*tot_site] = (double)System->CDown_2_D_RE[site]->Row[System->CDown_2_D_RE[site]->row_dim]/System->CDown_2_D_RE[site]->max_val;
   }
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      for (site = 0; site < tot_site; site++) {
         P[site + 13*tot_site]  = (double)System->SzL_LE[site]->Row[System->SzL_LE[site]->row_dim]/System->SzL_LE[site]->max_val;
         P[site + 14*tot_site] = (double)System->SpL_LE[site]->Row[System->SpL_LE[site]->row_dim]/System->SpL_LE[site]->max_val;
         P[site + 15*tot_site] = (double)System->SmL_LE[site]->Row[System->SmL_LE[site]->row_dim]/System->SmL_LE[site]->max_val;
         P[site + 16*tot_site] = (double)System->CUp_1_LE[site]->Row[System->CUp_1_LE[site]->row_dim]/System->CUp_1_LE[site]->max_val;
         P[site + 17*tot_site] = (double)System->CDown_1_LE[site]->Row[System->CDown_1_LE[site]->row_dim]/System->CDown_1_LE[site]->max_val;
         P[site + 18*tot_site] = (double)System->CUp_1_D_LE[site]->Row[System->CUp_1_D_LE[site]->row_dim]/System->CUp_1_D_LE[site]->max_val;
         P[site + 19*tot_site] = (double)System->CDown_1_D_LE[site]->Row[System->CDown_1_D_LE[site]->row_dim]/System->CDown_1_D_LE[site]->max_val;
         P[site + 20*tot_site] = (double)System->CUp_2_LE[site]->Row[System->CUp_2_LE[site]->row_dim]/System->CUp_2_LE[site]->max_val;
         P[site + 21*tot_site] = (double)System->CDown_2_LE[site]->Row[System->CDown_2_LE[site]->row_dim]/System->CDown_2_LE[site]->max_val;
         P[site + 22*tot_site] = (double)System->CUp_2_D_LE[site]->Row[System->CUp_2_D_LE[site]->row_dim]/System->CUp_2_D_LE[site]->max_val;
         P[site + 23*tot_site] = (double)System->CDown_2_D_LE[site]->Row[System->CDown_2_D_LE[site]->row_dim]/System->CDown_2_D_LE[site]->max_val;
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
