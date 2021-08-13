//
//  OUTPUT_ENERGY.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/29.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_ENERGY(MODEL_1DKLM_TVF *Model, EXACT_HAM_INFO *Ham_Info) {
   
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
      fprintf(file,"%1.1lf  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %+.15lf  %+.15lf  %e\n",
              (double)Model->spin_loc/2.0,
              Model->tot_parity,
              Model->t,
              Model->J,
              Model->I_xy,
              Model->I_z,
              Model->D_z,
              Model->h_x,
              Model->mu,
              Ham_Info->Value[num],
              Ham_Info->Value[num]/Model->tot_site,
              Ham_Info->Error[num]
              );
      fclose(file);
   }
}
