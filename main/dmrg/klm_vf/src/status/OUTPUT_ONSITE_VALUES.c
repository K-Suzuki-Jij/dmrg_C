//
//  OUTPUT_ONSITE_VALUES.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/07.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DKLM_VF *Model, DMRG_STATUS *Dmrg_Status) {
   
   int LL_site = Dmrg_Status->LL_site;
   int RR_site = Dmrg_Status->RR_site;
   int site;
   
   //Out put results
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/ExpectationValues", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/ExpectationValues/%s.txt", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now, Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_ONSITE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file,"###N=%d,Ne=%d,LocSpin=%1.1lf,Sz=%1.1lf,m=%d,BC=%s,Copy=%s,sweep=%d\n###t=%.1lf,J=%.4lf,I_xy=%.4lf,I_z=%.4lf,D_z=%.4lf,h_z=%.4lf,mu=%.4lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->spin_loc/2.0,
           (double)Model->tot_sz/2.0,
           Dmrg_Status->max_dim_system,
           Dmrg_Status->BC,
           Dmrg_Status->Enviro_Copy,
           Dmrg_Status->sweep_now,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_z,
           Model->mu
           );
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0) {
      int c_site = 0;
      for (site = 0; site <= LL_site + 1; site++) {
         fprintf(file,"%-3d  %+.15lf\n", c_site, Out[site]);
         c_site++;
      }
      for (site = Model->tot_site - 1; site >= LL_site + 2; site--) {
         fprintf(file,"%-3d  %+.15lf\n", c_site, Out[site]);
         c_site++;
      }
   }
   else {
      for (site = 0; site < Model->tot_site; site++) {
         fprintf(file,"%-3d  %+.15lf\n", site, Out[site]);
      }
   }
   fprintf(file,"\n");
   fclose(file);
   
}
