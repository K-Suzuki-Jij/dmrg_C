//
//  FREE_WHOLE_BASIS_Q1.c
//  1DKLM_TVF_EXACT
//
//  Created by Kohei Suzuki on 2019/07/29.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void FREE_WHOLE_BASIS_Q1(EXACT_WHOLE_BASIS_Q1 *W_Basis) {
   
   FREE_ARRAY_INT1(W_Basis->Dim);
   FREE_ARRAY_LINT2(W_Basis->Basis, 2);

}
