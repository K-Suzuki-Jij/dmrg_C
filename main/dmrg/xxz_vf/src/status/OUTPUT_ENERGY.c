//
//  OUTPUT_ENERGY.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/21.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_ENERGY(MODEL_1DXXZ_VF *Model, DMRG_STATUS *Dmrg_Status) {
   
   int LL_site = Dmrg_Status->LL_site;
   int RR_site = Dmrg_Status->RR_site;
   
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/AverageValues", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/AverageValues/energy.txt", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);

   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_ENERGY\n");
      printf("Can't open file\n");
      exit(1);
   }
   fprintf(file,"%1.1lf  %1.1lf  %.4lf  %.4lf  %.4lf  %.4lf  %+.15lf  %+.15lf  %d  %d\n",
           (double)Model->spin/2.0,
           (double)Model->tot_sz/2.0,
           Model->J_xy,
           Model->J_z,
           Model->D_z,
           Model->h_z,
           Dmrg_Status->gs_val,
           Dmrg_Status->gs_val/(LL_site + RR_site + 4),
           0,
           LL_site + RR_site + 3
           );
   fclose(file);
   
   
}
