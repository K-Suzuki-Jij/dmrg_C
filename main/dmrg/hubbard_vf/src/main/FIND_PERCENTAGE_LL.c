//
//  FIND_PERCENTAGE_LL.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DHUBBARD_VF *Model) {
   
   int tot_site  = Model->tot_site;
   int tot_num   = 30*tot_site;
   double *P     = GET_ARRAY_DOUBLE1(tot_num);
   int site;
   
   for (site = 0; site < tot_site; site++) {
      P[site + 0*tot_site] = (double)System->Ham[site]->Row[System->Ham[site]->row_dim]/System->Ham[site]->max_val;
      P[site + 1*tot_site] = (double)System->CUp_RE[site]->Row[System->CUp_RE[site]->row_dim]/System->CUp_RE[site]->max_val;
      P[site + 2*tot_site] = (double)System->CDown_RE[site]->Row[System->CDown_RE[site]->row_dim]/System->CDown_RE[site]->max_val;
      P[site + 3*tot_site] = (double)System->CUp_D_RE[site]->Row[System->CUp_D_RE[site]->row_dim]/System->CUp_D_RE[site]->max_val;
      P[site + 4*tot_site] = (double)System->CDown_D_RE[site]->Row[System->CDown_D_RE[site]->row_dim]/System->CDown_D_RE[site]->max_val;
      P[site + 5*tot_site] = (double)System->NC_Up_RE[site]->Row[System->NC_Up_RE[site]->row_dim]/System->NC_Up_RE[site]->max_val;
      P[site + 6*tot_site] = (double)System->NC_Down_RE[site]->Row[System->NC_Down_RE[site]->row_dim]/System->NC_Down_RE[site]->max_val;
   }
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      for (site = 0; site < tot_site; site++) {
         P[site + 7*tot_site]  = (double)System->NC_Up_LE[site]->Row[System->NC_Up_LE[site]->row_dim]/System->NC_Up_LE[site]->max_val;
         P[site + 8*tot_site] = (double)System->NC_Down_LE[site]->Row[System->NC_Down_LE[site]->row_dim]/System->NC_Down_LE[site]->max_val;
         P[site + 9*tot_site] = (double)System->CUp_LE[site]->Row[System->CUp_LE[site]->row_dim]/System->CUp_LE[site]->max_val;
         P[site + 10*tot_site] = (double)System->CDown_LE[site]->Row[System->CDown_LE[site]->row_dim]/System->CDown_LE[site]->max_val;
         P[site + 11*tot_site] = (double)System->CUp_D_LE[site]->Row[System->CUp_D_LE[site]->row_dim]/System->CUp_D_LE[site]->max_val;
         P[site + 12*tot_site] = (double)System->CDown_D_LE[site]->Row[System->CDown_D_LE[site]->row_dim]/System->CDown_D_LE[site]->max_val;
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
