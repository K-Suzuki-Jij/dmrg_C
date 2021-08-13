//
//  PRINT_MEM_STATUS.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/22.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void PRINT_MEM_STATUS(MODEL_1DTKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   int dim_onsite = Model->dim_onsite;
   int elem_num   = max_dim*max_dim*Dmrg_Param->sp_LL;
   
   double mem = 0;
   
   mem += pow(10,-9)*12*tot_site*(12*elem_num + 8*(max_dim + 1));
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      mem += pow(10,-9)*11*tot_site*(12*elem_num + 8*(max_dim + 1));
   }
   
   mem += pow(10,-9)*18*(12*dim_onsite*dim_onsite + 8*(dim_onsite + 1));
   
   mem += pow(10,-9)*4*4*tot_site*max_dim;
   
   mem += pow(10,-9)*8*(tot_site/2)*max_dim*dim_onsite;
   
   if (strcmp(Dmrg_Param->Enviro_Copy, "Yes") != 0) {
      mem = 2*mem;
   }
   
   printf("###N=%d,Ne=%d,%d,LocSpin=%1.1lf,Sz=%1.1lf,m=%d,BC=%s,Copy=%s,sweep=%d\n###t=%.1lf,J=%.4lf,I_xy=%.4lf,I_z=%.4lf,D_z=%.4lf,h_z=%.4lf,mu=%.4lf\n",
          Model->tot_site,
          Model->tot_ele_1,
          Model->tot_ele_2,
          (double)Model->spin_loc/2.0,
          (double)Model->tot_sz/2.0,
          max_dim,
          Model->BC,
          Dmrg_Param->Enviro_Copy,
          Dmrg_Param->sweep,
          Model->t,
          Model->J,
          Model->I_xy,
          Model->I_z,
          Model->D_z,
          Model->h_z,
          Model->mu
          );
   
   printf("Block_Mem=%lf[GB]\n",mem);
   
   mkdir("./result", 0777);
   FILE *file;
   if((file = fopen("./result/log.txt","a+")) == NULL){
      printf("Error in PRINT_STATUS\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file, "###N=%d,Ne=%d,%d,LocSpin=%1.1lf,Sz=%1.1lf,m=%d,BC=%s,Copy=%s,sweep=%d\n###t=%.1lf,J=%.4lf,I_xy=%.4lf,I_z=%.4lf,D_z=%.4lf,h_z=%.4lf,mu=%.4lf\n",
          Model->tot_site,
          Model->tot_ele_1,
          Model->tot_ele_2,
          (double)Model->spin_loc/2.0,
          (double)Model->tot_sz/2.0,
          max_dim,
          Model->BC,
          Dmrg_Param->Enviro_Copy,
          Dmrg_Param->sweep,
          Model->t,
          Model->J,
          Model->I_xy,
          Model->I_z,
          Model->D_z,
          Model->h_z,
          Model->mu
          );
   
   fprintf(file, "Block_Mem=%lf[GB]\n",mem);
   fclose(file);
   
   
}
