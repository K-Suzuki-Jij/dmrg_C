//
//  FIND_PERCENTAGE_LL.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

double FIND_PERCENTAGE_LL(BLOCK *System, MODEL_1DXXZ_VF *Model) {

   int tot_site  = Model->tot_site;
   int tot_num   = 12*tot_site;
   double *P     = GET_ARRAY_DOUBLE1(tot_num);
   int cf_length = Model->tot_site/2 - 1 - Model->cf_origin;
   int site;
   
   for (site = 0; site < tot_site; site++) {
      P[site + 0*tot_site] = (double)System->Ham[site]->Row[System->Ham[site]->row_dim]/System->Ham[site]->max_val;
      P[site + 1*tot_site] = (double)System->Sz_RE[site]->Row[System->Sz_RE[site]->row_dim]/System->Sz_RE[site]->max_val;
      P[site + 2*tot_site] = (double)System->Sp_RE[site]->Row[System->Sp_RE[site]->row_dim]/System->Sp_RE[site]->max_val;
      P[site + 3*tot_site] = (double)System->Sm_RE[site]->Row[System->Sm_RE[site]->row_dim]/System->Sm_RE[site]->max_val;
   }
   
   if (strcmp(Model->BC, "PBC") == 0) {
      for (site = 0; site < tot_site; site++) {
         P[site + 4*tot_site] = (double)System->Sz_LE[site]->Row[System->Sz_LE[site]->row_dim]/System->Sz_LE[site]->max_val;
         P[site + 5*tot_site] = (double)System->Sp_LE[site]->Row[System->Sp_LE[site]->row_dim]/System->Sp_LE[site]->max_val;
         P[site + 6*tot_site] = (double)System->Sm_LE[site]->Row[System->Sm_LE[site]->row_dim]/System->Sm_LE[site]->max_val;
      }
   }
   
   for (site = 0; site < tot_site/2; site++) {
      P[site + 7*tot_site] = (double)System->Sz[site]->Row[System->Sz[site]->row_dim]/System->Sz[site]->max_val;
      P[site + 8*tot_site] = (double)System->SzSz[site]->Row[System->SzSz[site]->row_dim]/System->SzSz[site]->max_val;
      P[site + 9*tot_site] = (double)System->Sx[site]->Row[System->Sx[site]->row_dim]/System->Sx[site]->max_val;
      P[site + 10*tot_site] = (double)System->SxSx[site]->Row[System->SxSx[site]->row_dim]/System->SxSx[site]->max_val;
   }
   
   for (site = 0; site < cf_length; site++) {
      P[site + 10*tot_site] = (double)System->Sz_CF[site]->Row[System->Sz_CF[site]->row_dim]/System->Sz_CF[site]->max_val;
      P[site + 11*tot_site] = (double)System->Sx_CF[site]->Row[System->Sx_CF[site]->row_dim]/System->Sx_CF[site]->max_val;
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
