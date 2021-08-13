//
//  FREE_HAM_BOX.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_HAM_BOX(HAM_BOX **Box, MODEL_1DHUBBARD_VF *Model) {
   
   int p_threads = Model->p_threads;
   
   int i;
   for (i = 0; i < p_threads; i++) {
      
      FREE_CRS1(Box[i]->Ham_On    );
      FREE_CRS1(Box[i]->CUp_On    );
      FREE_CRS1(Box[i]->CUp_D_On  );
      FREE_CRS1(Box[i]->CDown_On  );
      FREE_CRS1(Box[i]->CDown_D_On);
      FREE_CRS1(Box[i]->NC_Up_On  );
      FREE_CRS1(Box[i]->NC_Down_On);
      
      FREE_CRS1(Box[i]->Ham_System         );
      FREE_CRS1(Box[i]->CUp_RE_System      );
      FREE_CRS1(Box[i]->CDown_RE_System    );
      FREE_CRS1(Box[i]->CUp_D_RE_System    );
      FREE_CRS1(Box[i]->CDown_D_RE_System  );
      FREE_CRS1(Box[i]->NC_Up_RE_System    );
      FREE_CRS1(Box[i]->NC_Down_RE_System  );
      
      FREE_CRS1(Box[i]->Ham_Enviro         );
      FREE_CRS1(Box[i]->CUp_RE_Enviro      );
      FREE_CRS1(Box[i]->CDown_RE_Enviro    );
      FREE_CRS1(Box[i]->CUp_D_RE_Enviro    );
      FREE_CRS1(Box[i]->CDown_D_RE_Enviro  );
      FREE_CRS1(Box[i]->NC_Up_RE_Enviro    );
      FREE_CRS1(Box[i]->NC_Down_RE_Enviro  );
      
      if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
         FREE_CRS1(Box[i]->CUp_LE_System      );
         FREE_CRS1(Box[i]->CDown_LE_System    );
         FREE_CRS1(Box[i]->CUp_D_LE_System    );
         FREE_CRS1(Box[i]->CDown_D_LE_System  );
         FREE_CRS1(Box[i]->NC_Up_LE_System    );
         FREE_CRS1(Box[i]->NC_Down_LE_System  );
         
         FREE_CRS1(Box[i]->CUp_LE_Enviro      );
         FREE_CRS1(Box[i]->CDown_LE_Enviro    );
         FREE_CRS1(Box[i]->CUp_D_LE_Enviro    );
         FREE_CRS1(Box[i]->CDown_D_LE_Enviro  );
         FREE_CRS1(Box[i]->NC_Up_LE_Enviro    );
         FREE_CRS1(Box[i]->NC_Down_LE_Enviro  );
      }
      
      FREE_ARRAY_INT1(Box[i]->Ele_LL);
      FREE_ARRAY_INT1(Box[i]->Ele_RR);
      FREE_ARRAY_INT1(Box[i]->Ele_On);
      FREE_ARRAY_CHAR1(Box[i]->BC);
      
      free(Box[i]);
      
   }
   
   free(Box);
   
}
