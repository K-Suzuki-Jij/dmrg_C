//
//  OUTPUT_ONSITE_VALUES.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/29.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_ONSITE_VALUES(double *Out, char Name[], MODEL_1DKLM_TVF *Model) {
   
   //Out put results
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result/ExpectationValues");
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/ExpectationValues/%s.txt", Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_ONSITE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }
   
   fprintf(file, "###N=%d,Ne=%d,LocSpin=%1.1lf,P=%d,BC=%s\n###t=%.1lf,J=%.3lf,I_xy=%.3lf,I_z=%.3lf,D_z=%.3lf,h_x=%.3lf,mu=%.5lf\n",
           Model->tot_site,
           Model->tot_ele,
           (double)Model->spin_loc/2.0,
           Model->tot_parity,
           Model->BC,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_x,
           Model->mu
           );
   
   int site;
   
   for (site = 0; site < Model->tot_site; site++) {
      fprintf(file, "%-2d  %+.15lf\n", site, Out[site]);
   }
   
   fprintf(file,"\n");
   fclose(file);
   
}
