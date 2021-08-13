//
//  FREE_WHOLE_BASIS_Q3.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q3(EXACT_WHOLE_BASIS_Q3 *W_Basis, int max_sz) {
   
   FREE_ARRAY_INT3(W_Basis->Dim, 3, 2);
   FREE_ARRAY_LINT4(W_Basis->Basis, 3, 2, max_sz + 1);
   free(W_Basis);
   
}
