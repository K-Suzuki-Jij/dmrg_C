//
//  MAKE_ELEMENT_HAM_LLLR.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box) {
   
   int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_LL_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_System, 1.0, &elem_num);
   DMRG_MAKE_ELEM_LR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0*Box->SSD_LR, &elem_num);
   
   //Interaction LLLR
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->CUp_D_RE_System  , Box->CUp_On    , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->CUp_RE_System    , Box->CUp_D_On  , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->CDown_D_RE_System, Box->CDown_On  , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->CDown_RE_System  , Box->CDown_D_On, Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->NC_Up_RE_System  , Box->NC_Up_On  , Box->V*Box->SSD_LLLR , NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->NC_Down_RE_System, Box->NC_Down_On, Box->V*Box->SSD_LLLR , NULL       , &elem_num, NULL   , "No" );
   
   int num,LL,LR,inv;
   for (num = 0; num < elem_num; num++) {
      LL = A_Basis->LL_LLLRRRRL[num];
      LR = A_Basis->LR_LLLRRRRL[num];
      inv = Inv[LL][LR];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
