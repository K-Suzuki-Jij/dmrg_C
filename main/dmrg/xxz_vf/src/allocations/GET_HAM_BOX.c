//
//  GET_HAM_BOX.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/18.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DXXZ_VF *Model, int LL_site, int RR_site) {
   
   int p_threads  = Model->p_threads;
   int tot_site   = Model->tot_site;
   int i;
   
   HAM_BOX **Box = malloc(sizeof(HAM_BOX*)*p_threads);
   for (i = 0; i < p_threads; i++) {
      Box[i] = malloc(sizeof(HAM_BOX));
      Box[i]->Sp_On  = GET_COPY_CRS1(System->Sp_On );
      Box[i]->Sm_On  = GET_COPY_CRS1(System->Sm_On );
      Box[i]->Sz_On  = GET_COPY_CRS1(System->Sz_On );
      Box[i]->Ham_On = GET_COPY_CRS1(System->Ham_On);
      
      Box[i]->Ham_System   = GET_COPY_CRS1(System->Ham[LL_site]  );
      Box[i]->Sp_RE_System = GET_COPY_CRS1(System->Sp_RE[LL_site]);
      Box[i]->Sm_RE_System = GET_COPY_CRS1(System->Sm_RE[LL_site]);
      Box[i]->Sz_RE_System = GET_COPY_CRS1(System->Sz_RE[LL_site]);
      
      Box[i]->Ham_Enviro   = GET_COPY_CRS1(Enviro->Ham[RR_site]  );
      Box[i]->Sp_RE_Enviro = GET_COPY_CRS1(Enviro->Sp_RE[RR_site]);
      Box[i]->Sm_RE_Enviro = GET_COPY_CRS1(Enviro->Sm_RE[RR_site]);
      Box[i]->Sz_RE_Enviro = GET_COPY_CRS1(Enviro->Sz_RE[RR_site]);
      
      if (strcmp(Model->BC, "PBC") == 0) {
         Box[i]->Sp_LE_System = GET_COPY_CRS1(System->Sp_LE[LL_site]);
         Box[i]->Sm_LE_System = GET_COPY_CRS1(System->Sm_LE[LL_site]);
         Box[i]->Sz_LE_System = GET_COPY_CRS1(System->Sz_LE[LL_site]);
         
         Box[i]->Sp_LE_Enviro = GET_COPY_CRS1(Enviro->Sp_LE[RR_site]);
         Box[i]->Sm_LE_Enviro = GET_COPY_CRS1(Enviro->Sm_LE[RR_site]);
         Box[i]->Sz_LE_Enviro = GET_COPY_CRS1(Enviro->Sz_LE[RR_site]);
      }
      
      Box[i]->SSD_LR   = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR", "Onsite");
      Box[i]->SSD_RL   = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL", "Onsite");
      Box[i]->SSD_LLLR = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LL", "Intersite");
      Box[i]->SSD_LRRL = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR", "Intersite");
      Box[i]->SSD_RRRL = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL", "Intersite");
      
      Box[i]->J_xy = Model->J_xy;
      Box[i]->J_z  = Model->J_z ;
      Box[i]->D_z  = Model->D_z ;
      Box[i]->h_z  = Model->h_z ;

      Box[i]->BC   = GET_ARRAY_CHAR1(100);
      strcpy(Box[i]->BC, Model->BC);
      
   }
   
   return Box;
   
}
