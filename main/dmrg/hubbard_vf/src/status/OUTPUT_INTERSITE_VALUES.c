//
//  OUTPUT_INTERSITE_VALUES.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, char Name[], MODEL_1DHUBBARD_VF *Model, DMRG_STATUS *Dmrg_Status) {
   
   int LL_site = Dmrg_Status->LL_site;
   int RR_site = Dmrg_Status->RR_site;
   int origin = Model->cf_origin;
   int site,r;
   
   //Out put results
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/CorrelationFunctions", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/CorrelationFunctions/%s.txt", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now, Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_INTERSITE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file,"###N=%d,Ne=%d,Sz=%1.1lf,m=%d,BC=%s,Copy=%s,sweep=%d\n###t1=%.1lf,t2=%.4lf,U=%.4lf,V=%.4lf,h_z=%.4lf,mu=%.4lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->tot_sz/2.0,
           Dmrg_Status->max_dim_system,
           Dmrg_Status->BC,
           Dmrg_Status->Enviro_Copy,
           Dmrg_Status->sweep_now,
           Model->t1,
           Model->t2,
           Model->U,
           Model->V,
           Model->h_z,
           Model->mu
           );
   
   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0) {
      int c_site = 0;
      for (site = origin; site <= LL_site + 1; site++) {
         r = site - origin;
         fprintf(file,"%-3d  %-3d  %-3d  %+.15lf  %+.15lf\n",origin, site, r, Out[r], Out[r] - Onsite_Val[origin]*Onsite_Val[site]);
         c_site++;
      }
      for (site = Model->tot_site - 1; site >= LL_site + 2; site--) {
         r = site - origin;
         fprintf(file,"%-3d  %-3d  %-3d  %+.15lf  %+.15lf\n",origin, c_site + origin, c_site, Out[r], Out[r] - Onsite_Val[origin]*Onsite_Val[site]);
         c_site++;
      }
   }
   else {
      for (site = origin; site < Model->tot_site; site++) {
         r = site - origin;
         fprintf(file,"%-3d  %-3d  %-3d  %+.15lf  %+.15lf\n",origin, site, r, Out[r], Out[r] - Onsite_Val[origin]*Onsite_Val[site]);
      }
   }
   fprintf(file,"\n");
   fclose(file);
   
   
   
   
}
