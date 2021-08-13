//
//  FREE_HAM_BOX.c
//  1DTKLM_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/21.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DTKLM_VF *Model) {
   
   int p_threads  = Model->p_threads;
   
   int i;
   for (i = 0; i < p_threads; i++) {
      FREE_CRS1(Box[i]->Ham_On      );
      FREE_CRS1(Box[i]->SpL_On      );
      FREE_CRS1(Box[i]->SmL_On      );
      FREE_CRS1(Box[i]->SzL_On      );
      FREE_CRS1(Box[i]->CUp_1_On    );
      FREE_CRS1(Box[i]->CUp_1_D_On  );
      FREE_CRS1(Box[i]->CDown_1_On  );
      FREE_CRS1(Box[i]->CDown_1_D_On);
      FREE_CRS1(Box[i]->CUp_2_On    );
      FREE_CRS1(Box[i]->CUp_2_D_On  );
      FREE_CRS1(Box[i]->CDown_2_On  );
      FREE_CRS1(Box[i]->CDown_2_D_On);
      
      FREE_CRS1(Box[i]->Ham_System         );
      FREE_CRS1(Box[i]->SpL_RE_System      );
      FREE_CRS1(Box[i]->SmL_RE_System      );
      FREE_CRS1(Box[i]->SzL_RE_System      );
      FREE_CRS1(Box[i]->CUp_1_RE_System    );
      FREE_CRS1(Box[i]->CDown_1_RE_System  );
      FREE_CRS1(Box[i]->CUp_1_D_RE_System  );
      FREE_CRS1(Box[i]->CDown_1_D_RE_System);
      FREE_CRS1(Box[i]->CUp_2_RE_System    );
      FREE_CRS1(Box[i]->CDown_2_RE_System  );
      FREE_CRS1(Box[i]->CUp_2_D_RE_System  );
      FREE_CRS1(Box[i]->CDown_2_D_RE_System);

      FREE_CRS1(Box[i]->Ham_Enviro         );
      FREE_CRS1(Box[i]->SpL_RE_Enviro      );
      FREE_CRS1(Box[i]->SmL_RE_Enviro      );
      FREE_CRS1(Box[i]->SzL_RE_Enviro      );
      FREE_CRS1(Box[i]->CUp_1_RE_Enviro    );
      FREE_CRS1(Box[i]->CDown_1_RE_Enviro  );
      FREE_CRS1(Box[i]->CUp_1_D_RE_Enviro  );
      FREE_CRS1(Box[i]->CDown_1_D_RE_Enviro);
      FREE_CRS1(Box[i]->CUp_2_RE_Enviro    );
      FREE_CRS1(Box[i]->CDown_2_RE_Enviro  );
      FREE_CRS1(Box[i]->CUp_2_D_RE_Enviro  );
      FREE_CRS1(Box[i]->CDown_2_D_RE_Enviro);
      
      if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
         
         FREE_CRS1(Box[i]->SpL_LE_System      );
         FREE_CRS1(Box[i]->SmL_LE_System      );
         FREE_CRS1(Box[i]->SzL_LE_System      );
         FREE_CRS1(Box[i]->CUp_1_LE_System    );
         FREE_CRS1(Box[i]->CDown_1_LE_System  );
         FREE_CRS1(Box[i]->CUp_1_D_LE_System  );
         FREE_CRS1(Box[i]->CDown_1_D_LE_System);
         FREE_CRS1(Box[i]->CUp_2_LE_System    );
         FREE_CRS1(Box[i]->CDown_2_LE_System  );
         FREE_CRS1(Box[i]->CUp_2_D_LE_System  );
         FREE_CRS1(Box[i]->CDown_2_D_LE_System);
         
         FREE_CRS1(Box[i]->SpL_LE_Enviro      );
         FREE_CRS1(Box[i]->SmL_LE_Enviro      );
         FREE_CRS1(Box[i]->SzL_LE_Enviro      );
         FREE_CRS1(Box[i]->CUp_1_LE_Enviro    );
         FREE_CRS1(Box[i]->CDown_1_LE_Enviro  );
         FREE_CRS1(Box[i]->CUp_1_D_LE_Enviro  );
         FREE_CRS1(Box[i]->CDown_1_D_LE_Enviro);
         FREE_CRS1(Box[i]->CUp_2_LE_Enviro    );
         FREE_CRS1(Box[i]->CDown_2_LE_Enviro  );
         FREE_CRS1(Box[i]->CUp_2_D_LE_Enviro  );
         FREE_CRS1(Box[i]->CDown_2_D_LE_Enviro);
         
      }
      
      FREE_ARRAY_INT1(Box[i]->Ele_LL);
      FREE_ARRAY_INT1(Box[i]->Ele_RR);
      FREE_ARRAY_INT1(Box[i]->Ele_On);
      
      FREE_ARRAY_CHAR1(Box[i]->BC);

      free(Box[i]);
   }
   
   free(Box);
   
}
