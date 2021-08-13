//
//  GET_HAM_BOX.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

HAM_BOX **GET_HAM_BOX(BLOCK *System, BLOCK *Enviro, MODEL_1DHUBBARD_VF *Model, int LL_site, int RR_site) {
   
   int dim_LL     = System->Ham[LL_site]->row_dim;
   int dim_RR     = Enviro->Ham[RR_site]->row_dim;
   int dim_onsite = System->Ham_On->row_dim;
   int p_threads  = Model->p_threads;
   int tot_site   = Model->tot_site;
   int i;
   
   HAM_BOX **Box = malloc(sizeof(HAM_BOX*)*p_threads);
   for (i = 0; i < p_threads; i++) {
      Box[i] = malloc(sizeof(HAM_BOX));
      
      Box[i]->Ham_On       = GET_COPY_CRS1(System->Ham_On    );
      Box[i]->CUp_On       = GET_COPY_CRS1(System->CUp_On    );
      Box[i]->CUp_D_On     = GET_COPY_CRS1(System->CUp_D_On  );
      Box[i]->CDown_On     = GET_COPY_CRS1(System->CDown_On  );
      Box[i]->CDown_D_On   = GET_COPY_CRS1(System->CDown_D_On);
      Box[i]->NC_Up_On     = GET_COPY_CRS1(System->NC_Up_On  );
      Box[i]->NC_Down_On   = GET_COPY_CRS1(System->NC_Down_On);

      Box[i]->Ham_System          = GET_COPY_CRS1(System->Ham[LL_site]         );
      Box[i]->CUp_RE_System       = GET_COPY_CRS1(System->CUp_RE[LL_site]      );
      Box[i]->CDown_RE_System     = GET_COPY_CRS1(System->CDown_RE[LL_site]    );
      Box[i]->CUp_D_RE_System     = GET_COPY_CRS1(System->CUp_D_RE[LL_site]    );
      Box[i]->CDown_D_RE_System   = GET_COPY_CRS1(System->CDown_D_RE[LL_site]  );
      Box[i]->NC_Up_RE_System     = GET_COPY_CRS1(System->NC_Up_RE[LL_site]    );
      Box[i]->NC_Down_RE_System   = GET_COPY_CRS1(System->NC_Down_RE[LL_site]  );
      
      Box[i]->Ham_Enviro          = GET_COPY_CRS1(Enviro->Ham[RR_site]         );
      Box[i]->CUp_RE_Enviro       = GET_COPY_CRS1(Enviro->CUp_RE[RR_site]      );
      Box[i]->CDown_RE_Enviro     = GET_COPY_CRS1(Enviro->CDown_RE[RR_site]    );
      Box[i]->CUp_D_RE_Enviro     = GET_COPY_CRS1(Enviro->CUp_D_RE[RR_site]    );
      Box[i]->CDown_D_RE_Enviro   = GET_COPY_CRS1(Enviro->CDown_D_RE[RR_site]  );
      Box[i]->NC_Up_RE_Enviro     = GET_COPY_CRS1(System->NC_Up_RE[RR_site]    );
      Box[i]->NC_Down_RE_Enviro   = GET_COPY_CRS1(System->NC_Down_RE[RR_site]  );
      
      if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
         
         Box[i]->CUp_LE_System       = GET_COPY_CRS1(System->CUp_LE[LL_site]      );
         Box[i]->CDown_LE_System     = GET_COPY_CRS1(System->CDown_LE[LL_site]    );
         Box[i]->CUp_D_LE_System     = GET_COPY_CRS1(System->CUp_D_LE[LL_site]    );
         Box[i]->CDown_D_LE_System   = GET_COPY_CRS1(System->CDown_D_LE[LL_site]  );
         Box[i]->NC_Up_LE_System     = GET_COPY_CRS1(System->NC_Up_LE[LL_site]    );
         Box[i]->NC_Down_LE_System   = GET_COPY_CRS1(System->NC_Down_LE[LL_site]  );

         Box[i]->CUp_LE_Enviro       = GET_COPY_CRS1(Enviro->CUp_LE[RR_site]      );
         Box[i]->CDown_LE_Enviro     = GET_COPY_CRS1(Enviro->CDown_LE[RR_site]    );
         Box[i]->CUp_D_LE_Enviro     = GET_COPY_CRS1(Enviro->CUp_D_LE[RR_site]    );
         Box[i]->CDown_D_LE_Enviro   = GET_COPY_CRS1(Enviro->CDown_D_LE[RR_site]  );
         Box[i]->NC_Up_LE_Enviro     = GET_COPY_CRS1(System->NC_Up_LE[RR_site]    );
         Box[i]->NC_Down_LE_Enviro   = GET_COPY_CRS1(System->NC_Down_LE[RR_site]  );
      
      }
      
      Box[i]->Ele_LL = GET_ARRAY_INT1(dim_LL);
      Box[i]->Ele_RR = GET_ARRAY_INT1(dim_RR);
      Box[i]->Ele_On = GET_ARRAY_INT1(dim_onsite);
      
      COPY_INT1(System->Tot_Ele[LL_site], Box[i]->Ele_LL, dim_LL    , 1);
      COPY_INT1(Enviro->Tot_Ele[RR_site], Box[i]->Ele_RR, dim_RR    , 1);
      COPY_INT1(System->Tot_Ele[0]      , Box[i]->Ele_On, dim_onsite, 1);

      Box[i]->SSD_LR    = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR"  , "Onsite");
      Box[i]->SSD_RL    = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL"  , "Onsite");
      Box[i]->SSD_LLLR  = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LL"  , "Intersite");
      Box[i]->SSD_LRRL  = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "LR"  , "Intersite");
      Box[i]->SSD_RRRL  = DMRG_SSD_COEFF(LL_site, tot_site, Model->BC, "RL"  , "Intersite");

      Box[i]->t1  = Model->t1 ;
      Box[i]->t2  = Model->t2 ;
      Box[i]->U   = Model->U  ;
      Box[i]->V   = Model->V  ;
      Box[i]->h_z = Model->h_z;
      Box[i]->mu  = Model->mu ;
      
      Box[i]->BC   = GET_ARRAY_CHAR1(100);
      strcpy(Box[i]->BC, Model->BC);
      
   }
   
   return Box;
}
