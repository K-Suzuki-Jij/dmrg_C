//
//  OUTPUT_FOURIER_COMPONENTS.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_FOURIER_COMPONENTS(double *Val, char Name[], int start, int end, MODEL_1DTKLM_VF *Model, DMRG_STATUS *Dmrg_Status) {
   
   int site1,site2,r;
   int length = end - start;
   double temp_cos,temp_sin;
   double Wave[length];
   double Out[length];
   
      for (site1 = start; site1 < end; site1++) {
      temp_cos = 0;
      temp_sin = 0;
      r = site1 - start;
      Wave[r] = 2*M_PI*r/length;
      for (site2 = start; site2 < end; site2++) {
         temp_cos = temp_cos + Val[site2]*cos(Wave[r]*(site2 - start));
         temp_sin = temp_sin + Val[site2]*sin(Wave[r]*(site2 - start));
      }
      temp_cos = temp_cos*temp_cos;
      temp_sin = temp_sin*temp_sin;
      Out[r] = sqrt(temp_sin + temp_cos)/(end - start);
   }
   
   int LL_site = Dmrg_Status->LL_site;
   int RR_site = Dmrg_Status->RR_site;
   
   //Out put results
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/Fourier", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now);
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/[%d_1_1_%d]_%d/Fourier/%s.txt", LL_site + 1, RR_site + 1, Dmrg_Status->sweep_now, Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_ONSITE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file,"###N=%d,Ne=%d,%d,LocSpin=%1.1lf,Sz=%1.1lf,m=%d,BC=%s,Copy=%s,sweep=%d\n###t=%.1lf,J=%.4lf,I_xy=%.4lf,I_z=%.4lf,D_z=%.4lf,h_z=%.4lf,mu=%.4lf\n",
           Model->tot_site,
           Model->tot_ele_1,
           Model->tot_ele_2,
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
   
   for (site1 = start; site1 < end; site1++) {
      r = site1 - start;
      fprintf(file,"%-3d  %lf  %+.15lf  %d  %d\n", r, Wave[r], Out[r], start, end);
   }
   fprintf(file, "\n");
   fclose(file);
   
}

