//
//  main.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int main(void) {
   
   DMRG_PARAMETER *Dmrg_Param = malloc(sizeof(*Dmrg_Param));
   MODEL_1DXXZ_VF *Model      = malloc(sizeof(*Model));
   ////////////////////////////////////////////////
   strcpy(Dmrg_Param->Initial_Guess, "No" );
   strcpy(Dmrg_Param->Enviro_Copy  , "Yes");
   strcpy(Dmrg_Param->Diag_Method  , "Lanczos_Slow");
   strcpy(Dmrg_Param->Lan_Con      , "Normal");
   strcpy(Dmrg_Param->II_Type      , "CG");
   ////////////////////////////////////////////////
   Dmrg_Param->sweep = 3;
   Dmrg_Param->max_dim_system = 800;
   Dmrg_Param->sp_LL   = 0.6;
   Dmrg_Param->param_iter = 1;
   Dmrg_Param->dim_LLLRRRRL_limit = INT_MAX;
   ////////////////////////////////////////////////
   strcpy(Model->BC, "OBC");
   Model->p_threads = 4;
   Model->spin      = 1; //2*S
   Model->tot_sz    = 0; //2*tot_Sz
   Model->tot_site  = 20;
   Model->cf_origin = 0;
   ////////////////////////////////////////////////
   Model->J_xy = 1.0;
   Model->J_z  = 1.0;
   Model->D_z  = 0.0;
   Model->h_z  = 0.0;
   ////////////////////////////////////////////////
   Dmrg_Param->diag_acc      = pow(10,-14);
   Dmrg_Param->diag_min_step = 1;
   Dmrg_Param->diag_max_step = 800;
   ////////////////////////////////////////////////
   Dmrg_Param->cg_acc      = pow(10,-6);
   Dmrg_Param->cg_max_step = 1000;
   ////////////////////////////////////////////////
   Dmrg_Param->inv_iter_acc      = pow(10,-9);
   Dmrg_Param->inv_iter_diag_add = pow(10,-11);
   Dmrg_Param->inv_iter_max_step = 3;
   ////////////////////////////////////////////////
   
   for (Dmrg_Param->param_iter_now = 0; Dmrg_Param->param_iter_now < Dmrg_Param->param_iter; Dmrg_Param->param_iter_now++) {
      DMRG(Model, Dmrg_Param);
      Model->tot_sz = Model->tot_sz + 2;
   }
   
   free(Dmrg_Param);
   free(Model);
   
   return 0;
}
