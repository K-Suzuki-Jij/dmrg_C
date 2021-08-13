//
//  FREE_WHOLE_BASIS_Q2.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/30.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q2(EXACT_WHOLE_BASIS_Q2 *W_Basis) {
   
   FREE_ARRAY_INT2(W_Basis->Dim, 3);
   FREE_ARRAY_LINT3(W_Basis->Basis, 3, 2);
   free(W_Basis);
   
}
