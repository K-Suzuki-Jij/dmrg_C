//
//  OUTPUT_T_ERROR.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/14.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DTKLM_VF *Model) {
   
   int LL_site         = Dmrg_Status->LL_site;
   int system_size     = Dmrg_Status->LL_site + 2;
   int superblock_size = Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4;
   int dim_LLLR        = Dmrg_System->dim_LLLR;
   double ee,temp;
   int iter;
   char Out_Name[200];
   FILE *file;
   
   mkdir("./result", 0777);
   mkdir("./result/SystemInfo", 0777);
   
   //Entanglement Entropy
   sprintf(Out_Name, "./result/SystemInfo/EE_%d.txt", Dmrg_Status->sweep_now);
   
   ee = 0;
   for (iter = 0; iter < dim_LLLR; iter++) {
      temp = Dmrg_System->Val_DM_Dist[iter];
      if (temp > 0.0) {
         ee = ee + temp*log(temp);
      }
   }
   
   if((file = fopen(Out_Name,"a+")) == NULL){
      printf("Error in OUTPUT_T_ERROR\n");
      printf("Can't open file\n");
      exit(1);
   }
   
 if (LL_site == 0) {
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
   }
   fprintf(file,"%-3d  %-3d  %-3d  %+.15lf  %+.15lf\n",
           superblock_size,
           system_size,
           Dmrg_Status->tot_iter_now,
           log(superblock_size/M_PI*sin(M_PI*system_size/superblock_size)),
           -ee
           );
   if (LL_site + 5 == Model->tot_site || Dmrg_Status->tot_iter_now == Dmrg_Status->tot_iter) {
      fprintf(file,"\n");
   }
   fclose(file);
   
   //Trancation Error
   ee = 0;
   for (iter = 0; iter < dim_LLLR; iter++) {
      ee = ee + Dmrg_System->Val_DM_Dist[iter];;
   }
   
   sprintf(Out_Name, "./result/SystemInfo/T_Error_%d.txt", Dmrg_Status->sweep_now);
   file = fopen(Out_Name,"a+");
   if (LL_site == 0) {
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
   }
   fprintf(file, "%-3d  %-3d  %-3d  %+.15lf  %+.15lf\n", superblock_size, system_size, Dmrg_Status->tot_iter_now, ee, Dmrg_System->tr_error);
   if (LL_site + 5 == Model->tot_site || Dmrg_Status->tot_iter_now == Dmrg_Status->tot_iter) {
      fprintf(file,"\n");
   }
   fclose(file);
   
   //GS Energy
   sprintf(Out_Name, "./result/SystemInfo/Energy_%d.txt", Dmrg_Status->sweep_now);
   file = fopen(Out_Name, "a+");
   if (LL_site == 0) {
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
   }
   fprintf(file, "%-3d  %-3d  %-3d  %+.15lf  %+.15lf  %.1e\n",
           superblock_size, system_size, Dmrg_Status->tot_iter_now, Dmrg_Status->gs_val, Dmrg_Status->gs_val/superblock_size, Dmrg_Status->gs_error);
   if (LL_site + 5 == Model->tot_site || Dmrg_Status->tot_iter_now == Dmrg_Status->tot_iter) {
      fprintf(file,"\n");
   }
   fclose(file);
   
   
}
