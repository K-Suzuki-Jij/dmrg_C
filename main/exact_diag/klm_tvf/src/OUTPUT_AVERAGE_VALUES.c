//
//  OUTPUT_AVERAGE_VALUES.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/29.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_AVERAGE_VALUES(double *Out, char Name[], MODEL_1DKLM_TVF *Model) {
   
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result");
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/AverageValues");
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/AverageValues/avg_%s.txt", Name);
   
   FILE *file;
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_AVERAGE_VALUES\n");
      printf("Can't open file\n");
      exit(1);
   }

   int site;
   double val = 0;
   
   for (site = 0; site < Model->tot_site; site++) {
      val = val + Out[site];
   }
   
   fprintf(file,"%1.1lf  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %+.15lf  %+.15lf\n",
           (double)Model->spin_loc/2.0,
           Model->tot_parity,
           Model->t,
           Model->J,
           Model->I_xy,
           Model->I_z,
           Model->D_z,
           Model->h_x,
           Model->mu,
           val,
           val/Model->tot_site
           );
   
   fclose(file);
   
}
