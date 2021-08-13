//
//  GET_HAM_BOX.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DTKLM_VF *Model, int LL_site, int RR_site) {

   int dim_LL     = System->Ham[LL_site]->row_dim;
   int dim_RR     = Enviro->Ham[RR_site]->row_dim;
   int dim_onsite = System->Ham_On->row_dim;
   int p_threads  = Model->p_threads;
   int tot_site   = Model->tot_site;
   int i;
   
   HAM_BOX **Box = malloc(sizeof(HAM_BOX*)*p_threads);
   for (i = 0; i < p_threads; i++) {
      Box[i] = malloc(sizeof(HAM_BOX));
      Box[i]->Ham_On       = GET_COPY_CRS1(System->Ham_On      );
      Box[i]->SpL_On       = GET_COPY_CRS1(System->SpL_On      );
      Box[i]->SmL_On       = GET_COPY_CRS1(System->SmL_On      );
      Box[i]->SzL_On       = GET_COPY_CRS1(System->SzL_On      );
      Box[i]->CUp_1_On     = GET_COPY_CRS1(System->CUp_1_On    );
      Box[i]->CUp_1_D_On   = GET_COPY_CRS1(System->CUp_1_D_On  );
      Box[i]->CDown_1_On   = GET_COPY_CRS1(System->CDown_1_On  );
      Box[i]->CDown_1_D_On = GET_COPY_CRS1(System->CDown_1_D_On);
      Box[i]->CUp_2_On     = GET_COPY_CRS1(System->CUp_2_On    );
      Box[i]->CUp_2_D_On   = GET_COPY_CRS1(System->CUp_2_D_On  );
      Box[i]->CDown_2_On   = GET_COPY_CRS1(System->CDown_2_On  );
      Box[i]->CDown_2_D_On = GET_COPY_CRS1(System->CDown_2_D_On);
      
      Box[i]->Ham_System          = GET_COPY_CRS1(System->Ham[LL_site]         );
      Box[i]->SpL_RE_System       = GET_COPY_CRS1(System->SpL_RE[LL_site]      );
      Box[i]->SmL_RE_System       = GET_COPY_CRS1(System->SmL_RE[LL_site]      );
      Box[i]->SzL_RE_System       = GET_COPY_CRS1(System->SzL_RE[LL_site]      );
      Box[i]->CUp_1_RE_System     = GET_COPY_CRS1(System->CUp_1_RE[LL_site]    );
      Box[i]->CDown_1_RE_System   = GET_COPY_CRS1(System->CDown_1_RE[LL_site]  );
      Box[i]->CUp_1_D_RE_System   = GET_COPY_CRS1(System->CUp_1_D_RE[LL_site]  );
      Box[i]->CDown_1_D_RE_System = GET_COPY_CRS1(System->CDown_1_D_RE[LL_site]);
      Box[i]->CUp_2_RE_System     = GET_COPY_CRS1(System->CUp_2_RE[LL_site]    );
      Box[i]->CDown_2_RE_System   = GET_COPY_CRS1(System->CDown_2_RE[LL_site]  );
      Box[i]->CUp_2_D_RE_System   = GET_COPY_CRS1(System->CUp_2_D_RE[LL_site]  );
      Box[i]->CDown_2_D_RE_System = GET_COPY_CRS1(System->CDown_2_D_RE[LL_site]);
      
      Box[i]->Ham_Enviro          = GET_COPY_CRS1(Enviro->Ham[RR_site]         );
      Box[i]->SpL_RE_Enviro       = GET_COPY_CRS1(Enviro->SpL_RE[RR_site]      );
      Box[i]->SmL_RE_Enviro       = GET_COPY_CRS1(Enviro->SmL_RE[RR_site]      );
      Box[i]->SzL_RE_Enviro       = GET_COPY_CRS1(Enviro->SzL_RE[RR_site]      );
      Box[i]->CUp_1_RE_Enviro     = GET_COPY_CRS1(Enviro->CUp_1_RE[RR_site]    );
      Box[i]->CDown_1_RE_Enviro   = GET_COPY_CRS1(Enviro->CDown_1_RE[RR_site]  );
      Box[i]->CUp_1_D_RE_Enviro   = GET_COPY_CRS1(Enviro->CUp_1_D_RE[RR_site]  );
      Box[i]->CDown_1_D_RE_Enviro = GET_COPY_CRS1(Enviro->CDown_1_D_RE[RR_site]);
      Box[i]->CUp_2_RE_Enviro     = GET_COPY_CRS1(Enviro->CUp_2_RE[RR_site]    );
      Box[i]->CDown_2_RE_Enviro   = GET_COPY_CRS1(Enviro->CDown_2_RE[RR_site]  );
      Box[i]->CUp_2_D_RE_Enviro   = GET_COPY_CRS1(Enviro->CUp_2_D_RE[RR_site]  );
      Box[i]->CDown_2_D_RE_Enviro = GET_COPY_CRS1(Enviro->CDown_2_D_RE[RR_site]);
      
      if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
         
         Box[i]->SpL_LE_System       = GET_COPY_CRS1(System->SpL_LE[LL_site]      );
         Box[i]->SmL_LE_System       = GET_COPY_CRS1(System->SmL_LE[LL_site]      );
         Box[i]->SzL_LE_System       = GET_COPY_CRS1(System->SzL_LE[LL_site]      );
         Box[i]->CUp_1_LE_System     = GET_COPY_CRS1(System->CUp_1_LE[LL_site]    );
         Box[i]->CDown_1_LE_System   = GET_COPY_CRS1(System->CDown_1_LE[LL_site]  );
         Box[i]->CUp_1_D_LE_System   = GET_COPY_CRS1(System->CUp_1_D_LE[LL_site]  );
         Box[i]->CDown_1_D_LE_System = GET_COPY_CRS1(System->CDown_1_D_LE[LL_site]);
         Box[i]->CUp_2_LE_System     = GET_COPY_CRS1(System->CUp_2_LE[LL_site]    );
         Box[i]->CDown_2_LE_System   = GET_COPY_CRS1(System->CDown_2_LE[LL_site]  );
         Box[i]->CUp_2_D_LE_System   = GET_COPY_CRS1(System->CUp_2_D_LE[LL_site]  );
         Box[i]->CDown_2_D_LE_System = GET_COPY_CRS1(System->CDown_2_D_LE[LL_site]);
         
         Box[i]->SpL_LE_Enviro       = GET_COPY_CRS1(Enviro->SpL_LE[RR_site]      );
         Box[i]->SmL_LE_Enviro       = GET_COPY_CRS1(Enviro->SmL_LE[RR_site]      );
         Box[i]->SzL_LE_Enviro       = GET_COPY_CRS1(Enviro->SzL_LE[RR_site]      );
         Box[i]->CUp_1_LE_Enviro     = GET_COPY_CRS1(Enviro->CUp_1_LE[RR_site]    );
         Box[i]->CDown_1_LE_Enviro   = GET_COPY_CRS1(Enviro->CDown_1_LE[RR_site]  );
         Box[i]->CUp_1_D_LE_Enviro   = GET_COPY_CRS1(Enviro->CUp_1_D_LE[RR_site]  );
         Box[i]->CDown_1_D_LE_Enviro = GET_COPY_CRS1(Enviro->CDown_1_D_LE[RR_site]);
         Box[i]->CUp_2_LE_Enviro     = GET_COPY_CRS1(Enviro->CUp_2_LE[RR_site]    );
         Box[i]->CDown_2_LE_Enviro   = GET_COPY_CRS1(Enviro->CDown_2_LE[RR_site]  );
         Box[i]->CUp_2_D_LE_Enviro   = GET_COPY_CRS1(Enviro->CUp_2_D_LE[RR_site]  );
         Box[i]->CDown_2_D_LE_Enviro = GET_COPY_CRS1(Enviro->CDown_2_D_LE[RR_site]);
      }
      
      Box[i]->Ele_LL = GET_ARRAY_INT1(dim_LL);
      Box[i]->Ele_RR = GET_ARRAY_INT1(dim_RR);
      Box[i]->Ele_On = GET_ARRAY_INT1(dim_onsite);
      
      COPY_INT1(System->Tot_Ele[LL_site], Box[i]->Ele_LL, dim_LL    , 1);
      COPY_INT1(Enviro->Tot_Ele[RR_site], Box[i]->Ele_RR, dim_RR    , 1);
      COPY_INT1(System->Tot_Ele[0]      , Box[i]->Ele_On, dim_onsite, 1);

      Box[i]->SSD_LR   = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR", "Onsite");
      Box[i]->SSD_RL   = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL", "Onsite");
      Box[i]->SSD_LLLR = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LL", "Intersite");
      Box[i]->SSD_LRRL = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR", "Intersite");
      Box[i]->SSD_RRRL = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL", "Intersite");
      
      Box[i]->t    = Model->t   ;
      Box[i]->J    = Model->J   ;
      Box[i]->D_z  = Model->D_z ;
      Box[i]->I_xy = Model->I_xy;
      Box[i]->I_z  = Model->I_z ;
      Box[i]->h_z  = Model->h_z ;
      Box[i]->mu   = Model->mu  ;
      
      Box[i]->BC   = GET_ARRAY_CHAR1(100);
      strcpy(Box[i]->BC, Model->BC);
      
   }
   

   return Box;
}
