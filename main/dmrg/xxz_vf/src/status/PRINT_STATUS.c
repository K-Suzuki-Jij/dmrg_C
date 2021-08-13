//
//  PRINT_STATUS.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void PRINT_STATUS(DMRG_STATUS *Dmrg_Status, DMRG_TIME *Dmrg_Time) {
   
   Dmrg_Time->total = omp_get_wtime() - Dmrg_Time->total;
   
   if (strcmp(Dmrg_Status->BC, "OBC") == 0 || strcmp(Dmrg_Status->BC, "SSD") == 0) {
      printf("[%s]E/N=%.15lf(%.1e),N=%-3d(%-3d,1,1,%3d),dim=%-7d(%-4d,%-2d,%2d,%4d),Spin=%.1lf, Sz=%.1lf, T:%-6.1f(D:%-6.1f,II:%-6.1f,Ham:%-5.1f,Tr:%-4.1f),LL:%4.1lf%%,%3d/%3d, %2d/%d \n",
             Dmrg_Status->BC,
             Dmrg_Status->gs_val/(Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4),
             Dmrg_Status->gs_error,
             Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4,
             Dmrg_Status->LL_site + 1,
             Dmrg_Status->RR_site + 1,
             Dmrg_Status->dim_LLLRRRRL,
             Dmrg_Status->dim_LL,
             Dmrg_Status->dim_onsite,
             Dmrg_Status->dim_onsite,
             Dmrg_Status->dim_RR,
             Dmrg_Status->spin/2.0,
             Dmrg_Status->tot_sz/2.0,
             Dmrg_Time->total,
             Dmrg_Time->diag,
             Dmrg_Time->inv_iter,
             Dmrg_Time->make_ham,
             Dmrg_Time->trans_main + Dmrg_Time->trans_exp,
             Dmrg_Status->percent_LL*100,
             Dmrg_Status->tot_iter_now,
             Dmrg_Status->tot_iter,
             Dmrg_Status->param_iter_now + 1,
             Dmrg_Status->param_iter
             );
      
      mkdir("./result", 0777);
      FILE *file;
      if((file = fopen("./result/log.txt","a+")) == NULL){
         printf("Error in PRINT_STATUS\n");
         printf("Can't open file\n");
         exit(1);
      }
      
      fprintf(file, "[%s]E/N=%.15lf(%.1e),N=%-3d(%-3d,1,1,%3d),dim=%-7d(%-4d,%-2d,%2d,%4d),Spin=%.1lf, Sz=%.1lf, T:%-6.1f(D:%-6.1f,II:%-6.1f,Ham:%-5.1f,Tr:%-4.1f),LL:%4.1lf%%,%3d/%3d, %2d/%d \n",
              Dmrg_Status->BC,
              Dmrg_Status->gs_val/(Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4),
              Dmrg_Status->gs_error,
              Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4,
              Dmrg_Status->LL_site + 1,
              Dmrg_Status->RR_site + 1,
              Dmrg_Status->dim_LLLRRRRL,
              Dmrg_Status->dim_LL,
              Dmrg_Status->dim_onsite,
              Dmrg_Status->dim_onsite,
              Dmrg_Status->dim_RR,
              Dmrg_Status->spin/2.0,
              Dmrg_Status->tot_sz/2.0,
              Dmrg_Time->total,
              Dmrg_Time->diag,
              Dmrg_Time->inv_iter,
              Dmrg_Time->make_ham,
              Dmrg_Time->trans_main + Dmrg_Time->trans_exp,
              Dmrg_Status->percent_LL*100,
              Dmrg_Status->tot_iter_now,
              Dmrg_Status->tot_iter,
              Dmrg_Status->param_iter_now + 1,
              Dmrg_Status->param_iter
              );
      fclose(file);
   }
   
   else if (strcmp(Dmrg_Status->BC, "PBC") == 0) {
      printf("[%s]E/N=%.15lf(%.1e),N=%-3d(%-3d,1,%3d,1),dim=%-7d(%-4d,%-2d,%4d,%2d),Spin=%.1lf, Sz=%.1lf, T:%-6.1f(D:%-6.1f,II:%-6.1f,Ham:%-5.1f,Tr:%-4.1f),LL:%4.1lf%%,%3d/%3d, %2d/%d \n",
             Dmrg_Status->BC,
             Dmrg_Status->gs_val/(Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4),
             Dmrg_Status->gs_error,
             Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4,
             Dmrg_Status->LL_site + 1,
             Dmrg_Status->RR_site + 1,
             Dmrg_Status->dim_LLLRRRRL,
             Dmrg_Status->dim_LL,
             Dmrg_Status->dim_onsite,
             Dmrg_Status->dim_RR,
             Dmrg_Status->dim_onsite,
             Dmrg_Status->spin/2.0,
             Dmrg_Status->tot_sz/2.0,
             Dmrg_Time->total,
             Dmrg_Time->diag,
             Dmrg_Time->inv_iter,
             Dmrg_Time->make_ham,
             Dmrg_Time->trans_main + Dmrg_Time->trans_exp,
             Dmrg_Status->percent_LL*100,
             Dmrg_Status->tot_iter_now,
             Dmrg_Status->tot_iter,
             Dmrg_Status->param_iter_now + 1,
             Dmrg_Status->param_iter
             );
      
      mkdir("./result", 0777);
      FILE *file;
      if((file = fopen("./result/log.txt","a+")) == NULL){
         printf("Error in PRINT_STATUS\n");
         printf("Can't open file\n");
         exit(1);
      }
      
      fprintf(file, "[%s]E/N=%.15lf(%.1e),N=%-3d(%-3d,1,%3d,1),dim=%-7d(%-4d,%-2d,%4d,%2d),Spin=%.1lf, Sz=%.1lf, T:%-6.1f(D:%-6.1f,II:%-6.1f,Ham:%-5.1f,Tr:%-4.1f),LL:%4.1lf%%,%3d/%3d, %2d/%d \n",
              Dmrg_Status->BC,
              Dmrg_Status->gs_val/(Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4),
              Dmrg_Status->gs_error,
              Dmrg_Status->LL_site + Dmrg_Status->RR_site + 4,
              Dmrg_Status->LL_site + 1,
              Dmrg_Status->RR_site + 1,
              Dmrg_Status->dim_LLLRRRRL,
              Dmrg_Status->dim_LL,
              Dmrg_Status->dim_onsite,
              Dmrg_Status->dim_RR,
              Dmrg_Status->dim_onsite,
              Dmrg_Status->spin/2.0,
              Dmrg_Status->tot_sz/2.0,
              Dmrg_Time->total,
              Dmrg_Time->diag,
              Dmrg_Time->inv_iter,
              Dmrg_Time->make_ham,
              Dmrg_Time->trans_main + Dmrg_Time->trans_exp,
              Dmrg_Status->percent_LL*100,
              Dmrg_Status->tot_iter_now,
              Dmrg_Status->tot_iter,
              Dmrg_Status->param_iter_now + 1,
              Dmrg_Status->param_iter
              );
      fclose(file);
   }
   else {
      printf("Error in PRINT_STATUS\n");
      exit(1);
   }
   
   
}
