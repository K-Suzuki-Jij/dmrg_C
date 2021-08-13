//
//  FREE_HAM_BOX.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/21.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DXXZ_VF *Model) {
   
   int p_threads = Model->p_threads;
   
   int i;
   for (i = 0; i < p_threads; i++) {
      
      FREE_CRS1(Box[i]->Sp_On );
      FREE_CRS1(Box[i]->Sm_On );
      FREE_CRS1(Box[i]->Sz_On );
      FREE_CRS1(Box[i]->Ham_On);

      FREE_CRS1(Box[i]->Ham_System  );
      FREE_CRS1(Box[i]->Sp_RE_System);
      FREE_CRS1(Box[i]->Sm_RE_System);
      FREE_CRS1(Box[i]->Sz_RE_System);
      
      FREE_CRS1(Box[i]->Ham_Enviro  );
      FREE_CRS1(Box[i]->Sp_RE_Enviro);
      FREE_CRS1(Box[i]->Sm_RE_Enviro);
      FREE_CRS1(Box[i]->Sz_RE_Enviro);
      
      if (strcmp(Model->BC, "PBC") == 0) {
         FREE_CRS1(Box[i]->Sp_LE_System);
         FREE_CRS1(Box[i]->Sm_LE_System);
         FREE_CRS1(Box[i]->Sz_LE_System);
         
         FREE_CRS1(Box[i]->Sp_LE_Enviro);
         FREE_CRS1(Box[i]->Sm_LE_Enviro);
         FREE_CRS1(Box[i]->Sz_LE_Enviro);
      }
      
      FREE_ARRAY_CHAR1(Box[i]->BC);
      
      free(Box[i]);
   }
   
   free(Box);
   
}
