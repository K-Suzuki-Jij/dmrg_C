//
//  DMRG.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void DMRG(MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   Model->dim_onsite = 4;
   
   DMRG_STATUS *Dmrg_Status = GET_DMRG_STATUS(Model, Dmrg_Param);
   DMRG_TIME   *Dmrg_Time   = malloc(sizeof(DMRG_TIME));
   
   //PRINT_MEM_STATUS(Model, Dmrg_Param);
   
   BLOCK *System = GET_BLOCK_MATRIX(Model, Dmrg_Param);
   BLOCK *Enviro = System;

   if (strcmp(Dmrg_Param->Enviro_Copy, "Yes") != 0) {
      Enviro = GET_BLOCK_MATRIX(Model, Dmrg_Param);
   }
   
   //Infinite Algorithm
   for (Dmrg_Status->LL_site = 0; Dmrg_Status->LL_site < Model->tot_site/2 - 1; Dmrg_Status->LL_site++) {
      Dmrg_Status->RR_site = Dmrg_Status->LL_site;
      RENORMALIZE(System, System, Model, Dmrg_Param, Dmrg_Status, Dmrg_Time);
   }
   
   //Finite Algorithm
   for (Dmrg_Status->LL_site = Model->tot_site/2 - 1; Dmrg_Status->LL_site < Model->tot_site - 4; Dmrg_Status->LL_site++) {
      Dmrg_Status->RR_site = Model->tot_site - 4 - Dmrg_Status->LL_site;
      RENORMALIZE(System, System, Model, Dmrg_Param, Dmrg_Status, Dmrg_Time);
   }
   
   //SWEEP
   for (Dmrg_Status->sweep_now = 1; Dmrg_Status->sweep_now <= Dmrg_Param->sweep; Dmrg_Status->sweep_now++) {
      for (Dmrg_Status->LL_site = 0; Dmrg_Status->LL_site < Model->tot_site - 4; Dmrg_Status->LL_site++) {
         Dmrg_Status->RR_site = Model->tot_site - 4 - Dmrg_Status->LL_site;
         if (Dmrg_Status->LL_site == 0) {
            BLOCK *Temp = System;
            System = Enviro;
            Enviro = Temp;
         }
         RENORMALIZE(System, Enviro, Model, Dmrg_Param, Dmrg_Status, Dmrg_Time);
         if (Dmrg_Status->LL_site == Dmrg_Status->RR_site && Dmrg_Status->sweep_now == Dmrg_Param->sweep) {
            break;
         }
      }
   }
   
   free(Dmrg_Status);
   free(Dmrg_Time);
   
   FREE_BLOCK_MATRIX(System, Model, Dmrg_Param);
   if (strcmp(Dmrg_Param->Enviro_Copy, "Yes") != 0) {
      FREE_BLOCK_MATRIX(Enviro, Model, Dmrg_Param);
   }
   
}
