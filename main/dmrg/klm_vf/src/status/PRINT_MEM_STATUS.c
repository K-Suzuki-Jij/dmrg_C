//
//  PRINT_MEM_STATUS.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void PRINT_MEM_STATUS(MODEL_1DKLM_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   int dim_onsite = Model->dim_onsite;
   int elem_num   = max_dim*max_dim*Dmrg_Param->sp_LL;
   
   double mem = 0;
   
   mem += pow(10,-9)*8*tot_site*(12*elem_num + 8*(max_dim + 1));
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      mem += pow(10,-9)*7*tot_site*(12*elem_num + 8*(max_dim + 1));
   }
   
   mem += pow(10,-9)*2*tot_site*max_dim*4;
   
   mem += pow(10,-9)*(8 + 5)*(12*dim_onsite*dim_onsite + 8*(dim_onsite + 1));
   
   mem += pow(10,-9)*4*(tot_site/2)*max_dim*dim_onsite;

   if (strcmp(Dmrg_Param->Enviro_Copy, "Yes") != 0) {
      mem = 2*mem;
   }
   
   printf("t=%lf,J=%lf,D_z=%lf,I_xy=%lf,I_z=%lf,h_z=%lf,tot_site=%d,tot_ele=%d,local_spin=%.1lf,tot_sz=%.1lf,BC=%s,Enviro_Copy=%s\n",
          Model->t,
          Model->J,
          Model->D_z,
          Model->I_xy,
          Model->I_z,
          Model->h_z,
          Model->tot_site,
          Model->tot_ele,
          Model->spin_loc/2.0,
          Model->tot_sz/2.0,
          Model->BC,
          Dmrg_Param->Enviro_Copy
          );
   
   printf("Block_Mem=%lf[GB]\n",mem);
   
   mkdir("./result", 0777);
   FILE *file;
   if((file = fopen("./result/log.txt","a+")) == NULL){
      printf("Error in PRINT_STATUS\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file, "t=%lf,J=%lf,D_z=%lf,I_xy=%lf,I_z=%lf,h_z=%lf,tot_site=%d,tot_ele=%d,local_spin=%.1lf,tot_sz=%.1lf,BC=%s,Enviro_Copy=%s\n",
          Model->t,
          Model->J,
          Model->D_z,
          Model->I_xy,
          Model->I_z,
          Model->h_z,
          Model->tot_site,
          Model->tot_ele,
          Model->spin_loc/2.0,
          Model->tot_sz/2.0,
          Model->BC,
          Dmrg_Param->Enviro_Copy
          );
   
   fprintf(file, "Block_Mem=%lf[GB]\n",mem);
   fclose(file);
   
   
}
