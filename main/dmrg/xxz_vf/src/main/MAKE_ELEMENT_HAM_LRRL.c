//
//  MAKE_ELEMENT_HAM_LRRL.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/27.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box) {

    int elem_num = 0;
   
   //Interaction LRRL
   DMRG_MAKE_ELEM_LRRL_LRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_On, Box->Sm_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LRRL_LRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_On, Box->Sp_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LRRL_LRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_On, Box->Sz_On, Box->J_z     , NULL, &elem_num, NULL, "No");
   
   int num,LR,RL,inv;
   for (num = 0; num < elem_num; num++) {
      LR = A_Basis->LR_LLLRRRRL[num];
      RL = A_Basis->RL_LLLRRRRL[num];
      inv = Inv[LR][RL];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
