//
//  FREE_T_MATRIX.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_T_MATRIX(BLOCK *System, BLOCK *Enviro, int tot_site, int sweep_now, int sweep, char Enviro_Copy[100]) {
   
   int LL_site;
   for (LL_site = 0; LL_site <= tot_site/2 - 3; LL_site++) {
      FREE_CRS1(Enviro->TM_D[LL_site]);
      FREE_CCS1(Enviro->TM[LL_site]);
      FREE_ARRAY_SINT1(Enviro->Basis_LL_LLLR[LL_site]);
      FREE_ARRAY_SINT1(Enviro->Basis_LR_LLLR[LL_site]);
   }
   
   if (sweep_now == sweep && strcmp(Enviro_Copy, "Yes") != 0) {
      for (LL_site = 0; LL_site <= tot_site/2 - 3; LL_site++) {
         FREE_CRS1(System->TM_D[LL_site]);
         FREE_CCS1(System->TM[LL_site]);
         FREE_ARRAY_SINT1(System->Basis_LL_LLLR[LL_site]);
         FREE_ARRAY_SINT1(System->Basis_LR_LLLR[LL_site]);
      }
   }
   
}
