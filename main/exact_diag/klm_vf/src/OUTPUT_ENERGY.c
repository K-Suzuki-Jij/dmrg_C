//
//  OUTPUT_ENERGY.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/11.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_ENERGY(MODEL_1DKLM_VF *Model, EXACT_HAM_INFO *Ham_Info) {
   
   char Out_Name[200];
   mkdir("./result",0777);
   sprintf(Out_Name,"./result");
   mkdir(Out_Name,0777);
   sprintf(Out_Name,"./result/AverageValues");
   mkdir(Out_Name,0777);
   
   int num;
   
   for (num = 0; num < Ham_Info->diag_num; num++) {
      sprintf(Out_Name,"./result/AverageValues/energy%d.txt", num);
      
      FILE *file;
      if((file = fopen(Out_Name,"a+")) == NULL){
         printf("Error in OUTPUT_ENERGY\n");
         printf("Can't open file\n");
         exit(1);
      }
      fprintf(file,"%1.1lf  %1.1lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %+.15lf  %+.15lf  %e\n",
              (double)Model->spin_loc/2.0,
              (double)Model->tot_sz/2,
              Model->t,
              Model->J,
              Model->I_xy,
              Model->I_z,
              Model->D_z,
              Model->h_z,
              Model->h_z,
              Model->mu,
              Ham_Info->Value[num],
              Ham_Info->Value[num]/Model->tot_site,
              Ham_Info->Error[num]
              );
      fclose(file);
   }
}
