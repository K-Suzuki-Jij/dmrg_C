//
//  GET_DMRG_STATUS.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

DMRG_STATUS *GET_DMRG_STATUS(MODEL_1DKLM_TVF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   if (Model->dim_onsite >= Dmrg_Param->max_dim_system) {
      printf("Error in GET_DMRG_STATUS\n");
      printf("Need more max_dim = %d >= %d\n", Dmrg_Param->max_dim_system, Model->dim_onsite);
      exit(1);
   }
   
   DMRG_STATUS *Dmrg_Status    = malloc(sizeof(DMRG_STATUS));
   Dmrg_Status->param_iter     = Dmrg_Param->param_iter;
   Dmrg_Status->param_iter_now = Dmrg_Param->param_iter_now;
   Dmrg_Status->sweep          = Dmrg_Param->sweep;
   Dmrg_Status->sweep_now      = 0;
   Dmrg_Status->tot_iter_now   = 0;
   Dmrg_Status->tot_iter       = Model->tot_site/2 - 1 + Dmrg_Param->sweep*(Model->tot_site - 4);
   Dmrg_Status->dim_onsite     = Model->dim_onsite;
   Dmrg_Status->spin           = Model->spin_loc;
   Dmrg_Status->tot_parity     = Model->tot_parity;
   Dmrg_Status->tot_ele        = Model->tot_ele;
   Dmrg_Status->max_dim_system = Dmrg_Param->max_dim_system;
   
   strcpy(Dmrg_Status->BC, Model->BC);
   strcpy(Dmrg_Status->Enviro_Copy, Dmrg_Param->Enviro_Copy);
   
   return Dmrg_Status;
   
}
