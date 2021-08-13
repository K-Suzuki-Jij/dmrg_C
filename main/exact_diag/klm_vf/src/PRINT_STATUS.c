//
//  PRINT_STATUS.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void PRINT_STATUS(MODEL_1DKLM_VF *Model, EXACT_TIME *Time, EXACT_HAM_INFO *Ham_Info, EXACT_PARAMETER *Param) {
   
   printf("###N=%d,Ne=%d,LocSpin=%1.1lf,Sz=%1.1lf,BC=%s,Dim=%d,t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf,%d/%d\n",
          Model->tot_site,
          Model->tot_ele,
          (double)Model->spin_loc/2.0,
          (double)Model->tot_sz/2,
          Model->BC,
          Ham_Info->dim,
          Model->t,
          Model->J,
          Model->I_xy,
          Model->I_z,
          Model->D_z,
          Model->h_z,
          Model->mu,
          Param->param_iter_now + 1,
          Param->param_iter
          );
   
   printf("Mem: Basis=%lf, Ham=%lf [GB]\n",
          (double)Ham_Info->dim*8*pow(10,-9),
          Ham_Info->mem_ham
          );
   
   if (Param->param_iter_now >= 1) {
      Time->make_basis = 0;
      Time->find_dim = 0;
   }
   
   printf("Time=%.1lf: Diag=%.1lf, inv=%.1lf, ham=%.1lf, dim=%.1lf, basis=%.1lf, exp=%.1lf, cf=%.1lf, sc=%.1lf\n",
          Time->total,
          Time->diag,
          Time->inv_iter,
          Time->make_ham,
          Time->find_dim,
          Time->make_basis,
          Time->exp_values,
          Time->cf,
          Time->sc
          );
   
}
