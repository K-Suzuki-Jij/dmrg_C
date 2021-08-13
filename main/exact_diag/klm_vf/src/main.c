//
//  main.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int main(int argc, const char * argv[]) {
   
   MODEL_1DKLM_VF  *Model = malloc(sizeof(*Model));
   EXACT_PARAMETER *Param = malloc(sizeof(*Param));
   ////////////////////////////////////////////////
   strcpy(Model->BC, "PBC");
   Model->p_threads = 4;
   Model->tot_site  = 8;
   Model->tot_ele   = 8;
   Model->tot_sz    = 0; //2*tot_Sz
   Model->spin_loc  = 1; //2*S
   Model->cf_origin = 0;
   ////////////////////////////////////////////////
   Param->param_iter = 1;
   Model->t    = -1.0;
   Model->J    = +1.5;
   Model->D_z  = -1.0;
   Model->I_xy = -0.0;
   Model->I_z  = -1.0;
   Model->h_z  = +0.0;
   Model->mu   = +0.0;
   ////////////////////////////////////////////////
   strcpy(Param->Diag_Method, "Lanczos_Slow");
   Param->diag_acc      = pow(10,-14);
   Param->diag_min_step = 1;
   Param->diag_max_step = 800;
   Param->diag_num      = 1;
   ////////////////////////////////////////////////
   Param->cg_acc      = pow(10,-7);
   Param->cg_max_step = 1000;
   ////////////////////////////////////////////////
   Param->inv_iter_acc      = pow(10,-9);
   Param->inv_iter_diag_add = pow(10,-11);
   Param->inv_iter_max_step = 3;
   ////////////////////////////////////////////////
   Param->est_max_row_elem_num = 10000;
   ////////////////////////////////////////////////

   Model->dim_charge      = 4;
   Model->dim_lspin       = (Model->spin_loc + 1);
   Model->dim_onsite      = Model->dim_lspin*Model->dim_charge;
   Model->dim_csl_onsite  = Model->dim_lspin*Model->dim_lspin*2;
   Model->dim_ccsl_onsite = Model->dim_lspin*Model->dim_lspin;
   
   EXACT_BASIS_INFO *Basis_Info = malloc(sizeof(*Basis_Info));
   EXACT_HAM_INFO   *Ham_Info   = malloc(sizeof(*Ham_Info));
   EXACT_TIME       *Time       = malloc(sizeof(*Time));
   
   Basis_Info->dim   = FIND_DIM(Model, Time);
   Basis_Info->Basis = GET_BASIS(Model, Time, Basis_Info->dim);
   Ham_Info->Vector  = GET_ARRAY_DOUBLE2(Param->diag_num, Basis_Info->dim);
   Ham_Info->Value   = GET_ARRAY_DOUBLE1(Param->diag_num);
   Ham_Info->Error   = GET_ARRAY_DOUBLE1(Param->diag_num);
   
   for (Param->param_iter_now = 0; Param->param_iter_now < Param->param_iter; Param->param_iter_now++) {
      EXACT_DIAGONALIZATION(Model, Param, Basis_Info, Ham_Info, Time);
   }
   
   
   
   return 0;
}

