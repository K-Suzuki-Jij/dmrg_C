//
//  FREE_HAM_BOX.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/26.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_HAM_BOX(HAM_BOX *Box, MODEL_1DKLM_TVF *Model) {
   
   FREE_CRS1(Box->Ham_On   );
   FREE_CRS1(Box->SpL_On   );
   FREE_CRS1(Box->SmL_On   );
   FREE_CRS1(Box->SzL_On   );
   FREE_CRS1(Box->Even_On  );
   FREE_CRS1(Box->Even_D_On);
   FREE_CRS1(Box->Odd_On   );
   FREE_CRS1(Box->Odd_D_On );
   FREE_CRS1(Box->Zero_On  );
   FREE_CRS1(Box->SxL_On   );
   FREE_CRS1(Box->SxC_On   );
   FREE_CRS1(Box->SzC_On   );
   FREE_CRS1(Box->SzLSzL_On);
   FREE_CRS1(Box->SxLSxL_On);
   FREE_CRS1(Box->SzCSzC_On);
   FREE_CRS1(Box->NC_On    );
   FREE_CRS1(Box->NCNC_On  );
   FREE_CRS2(Box->CCSL_On, Model->dim_ccsl_onsite);
   FREE_CRS2(Box->CSL_On , Model->dim_csl_onsite );

   free(Box);

}
