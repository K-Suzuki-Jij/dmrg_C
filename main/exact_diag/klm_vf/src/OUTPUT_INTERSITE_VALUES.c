//
//  OUTPUT_INTERSITE_VALUES.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/11.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_INTERSITE_VALUES(double *Out, double *Onsite_Val, int origin, int end, char Name[], MODEL_1DKLM_VF *Model) {
   
   //Out put results
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/CorrelationFunctions");
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/CorrelationFunctions/%s.txt", Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_INTERSITE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file, "###N=%d,Ne=%d,LocSpin=%1.1lf,Sz=%1.1lf,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->spin_loc/2.0,
           (double)Model->tot_sz/2,
           Model->BC,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_z,
           Model->mu
           );
   
   int site,r;
   
   for (site = origin; site < end; site++) {
      r = site - origin;
      fprintf(file, "%-2d  %-2d  %-2d  %+.15lf  %+.15lf\n", origin, site, r, Out[r], Out[r] - Onsite_Val[origin]*Onsite_Val[site]);
   }
   
   fprintf(file,"\n");
   fclose(file);
   
}
