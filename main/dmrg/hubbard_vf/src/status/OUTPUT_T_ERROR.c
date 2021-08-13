//
//  OUTPUT_T_ERROR.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void OUTPUT_T_ERROR(DMRG_SYSTEM_INFO *Dmrg_System, DMRG_STATUS *Dmrg_Status, MODEL_1DHUBBARD_VF *Model) {
   
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
   }
   fprintf(file, "%-3d  %-3d  %-3d  %+.15lf  %+.15lf  %.1e\n",
           superblock_size, system_size, Dmrg_Status->tot_iter_now, Dmrg_Status->gs_val, Dmrg_Status->gs_val/superblock_size, Dmrg_Status->gs_error);
   if (LL_site + 5 == Model->tot_site || Dmrg_Status->tot_iter_now == Dmrg_Status->tot_iter) {
      fprintf(file,"\n");
   }
   fclose(file);
   
}
