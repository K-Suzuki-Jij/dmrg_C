//
//  FIND_SITE_STATE.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/09.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int FIND_SITE_STATE(long basis, int site, int dim_onsite) {
   
   int i;
   
   for (i = 0; i < site; i++) {
      basis = basis/(long)dim_onsite;
   }
   
   basis = basis%dim_onsite;
   
   return (int)basis;
   
}

