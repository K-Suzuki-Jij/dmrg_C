//
//  FREE_WHOLE_BASIS_Q2.c
//  1DKLM_VF_EXACT
//
//  Created by Kohei Suzuki on 2019/08/12.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2(EXACT_WHOLE_BASIS_Q2 *W_Basis, int max_sz) {
   
   FREE_ARRAY_INT2(W_Basis->Dim, 2);
   FREE_ARRAY_LINT3(W_Basis->Basis, 2, max_sz + 1);
   free(W_Basis);
   
}
