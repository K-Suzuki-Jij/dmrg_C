//
//  main.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int main(int argc, const char * argv[]) {
   
   DMRG_PARAMETER *Dmrg_Param = malloc(sizeof(*Dmrg_Param));
   MODEL_1DKLM_TVF *Model     = malloc(sizeof(*Model));
   ////////////////////////////////////////////////
   strcpy(Dmrg_Param->Initial_Guess, "No" );
   strcpy(Dmrg_Param->Enviro_Copy  , "Yes");
   strcpy(Dmrg_Param->Diag_Method  , "Lanczos_Slow");
   strcpy(Dmrg_Param->Lan_Con      , "Normal");
   strcpy(Dmrg_Param->II_Type      , "CG");
   ////////////////////////////////////////////////
   Dmrg_Param->sweep = 3;
   Dmrg_Param->dim_LLLRRRRL_limit = INT_MAX;
   Dmrg_Param->max_dim_system = 100;
   Dmrg_Param->sp_LL = 0.4;
   ////////////////////////////////////////////////
   strcpy(Model->BC, "OBC");
   Model->p_threads  = 4;
   Model->tot_site   = 20;
   Model->tot_ele    = 10;
   Model->tot_parity = 0; //0 or 1
   Model->spin_loc   = 2; //2*S
   Model->cf_origin  = 0;
   ////////////////////////////////////////////////
   Dmrg_Param->param_iter = 1;
   Model->t = -1.0;
   double J_Array[]    = {+2.0};
   double D_z_Array[]  = {-1.0};
   double I_xy_Array[] = {-0.3};
   double I_z_Array[]  = {-0.6};
   double h_xc_Array[] = {-0.5};
   double h_xl_Array[] = {-0.6};
   double mu_Array[]   = {+0.0};
   ////////////////////////////////////////////////
   Dmrg_Param->diag_acc      = pow(10,-14);
   Dmrg_Param->diag_min_step = 1;
   Dmrg_Param->diag_max_step = 800;
   ////////////////////////////////////////////////
   Dmrg_Param->cg_acc      = pow(10,-7);
   Dmrg_Param->cg_max_step = 1000;
   ////////////////////////////////////////////////
   Dmrg_Param->inv_iter_acc      = pow(10,-9);
   Dmrg_Param->inv_iter_diag_add = pow(10,-11);
   Dmrg_Param->inv_iter_max_step = 3;
   ////////////////////////////////////////////////

   for (Dmrg_Param->param_iter_now = 0; Dmrg_Param->param_iter_now < Dmrg_Param->param_iter; Dmrg_Param->param_iter_now++) {
      Model->J    = J_Array[Dmrg_Param->param_iter_now];
      Model->D_z  = D_z_Array[Dmrg_Param->param_iter_now];
      Model->I_xy = I_xy_Array[Dmrg_Param->param_iter_now];
      Model->I_z  = I_z_Array[Dmrg_Param->param_iter_now];
      Model->h_xc = h_xc_Array[Dmrg_Param->param_iter_now];
      Model->h_xl = h_xl_Array[Dmrg_Param->param_iter_now];
      Model->mu   = mu_Array[Dmrg_Param->param_iter_now];
      DMRG(Model, Dmrg_Param);
   }
   
   free(Model);
   free(Dmrg_Param);
   
   return 0;
}
