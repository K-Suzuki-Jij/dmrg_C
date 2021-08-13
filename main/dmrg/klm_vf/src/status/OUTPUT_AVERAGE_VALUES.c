//
//  OUTPUT_AVERAGE_VALUES.c
//  1DKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/07.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_AVERAGE_VALUES(double *Val, char Name[], int start, int end, MODEL_1DKLM_VF *Model, DMRG_STATUS *Dmrg_Status) {
 
   int site;
   double val = 0;
   
   for (site = start; site < end; site++) {
      val = val + Val[site];
   }
   
   int LL_site = Dmrg_Status->LL_site;
   int RR_site = Dmrg_Status->RR_site;
   
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/AverageValues", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/AverageValues/%s.txt", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now, Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_ENERGY\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file,"%1.1lf  %1.1lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %+.15lf  %d  %d\n",
           (double)Model->spin_loc/2.0,
           (double)Model->tot_sz/2.0,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_z,
           Model->mu,
           val/(end - start),
           start,
           end
           );
   fclose(file);
   
   
}
