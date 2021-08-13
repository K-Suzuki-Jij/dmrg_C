//
//  MAKE_ELEMENT_HAM_LLLR.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/18.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box) {

    int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_LL_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_System, 1.0, &elem_num);
   DMRG_MAKE_ELEM_LR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0, &elem_num);
   
   //Interaction LLLR
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Sp_RE_System, Box->Sm_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Sm_RE_System, Box->Sp_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Sz_RE_System, Box->Sz_On, Box->J_z , NULL, &elem_num, NULL, "No");
   
   int num,LL,LR,inv;
   for (num = 0; num < elem_num; num++) {
      LL = A_Basis->LL_LLLRRRRL[num];
      LR = A_Basis->LR_LLLRRRRL[num];
      inv = Inv[LL][LR];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
