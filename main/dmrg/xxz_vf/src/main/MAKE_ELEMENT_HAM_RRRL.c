//
//  MAKE_ELEMENT_HAM_RRRL.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/10/27.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_RRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box) {

    int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_RR_RRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_Enviro, 1.0, &elem_num);
   DMRG_MAKE_ELEM_RL_RRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0, &elem_num);
   
   //Interaction RRRL
   DMRG_MAKE_ELEM_RRRL_RRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_RE_Enviro, Box->Sm_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_RRRL_RRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_RE_Enviro, Box->Sp_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_RRRL_RRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_RE_Enviro, Box->Sz_On, Box->J_z , NULL, &elem_num, NULL, "No");
   
   int num,RR,RL,inv;
   for (num = 0; num < elem_num; num++) {
      RR = A_Basis->RR_LLLRRRRL[num];
      RL = A_Basis->RL_LLLRRRRL[num];
      inv = Inv[RR][RL];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
