//
//  MAKE_ELEMENT_HAM_LLLR.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LLLR(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int **Inv, HAM_BOX *Box) {
   
   int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_LL_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_System, 1.0, &elem_num);
   DMRG_MAKE_ELEM_LR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0*Box->SSD_LR, &elem_num);
   
   //Interaction LLLR
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Even_D_RE_System , Box->Even_On  , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Even_RE_System   , Box->Even_D_On, Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Odd_D_RE_System  , Box->Odd_On   , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->Odd_RE_System    , Box->Odd_D_On , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->SzL_RE_System    , Box->SzL_On   , Box->I_z*Box->SSD_LLLR     , NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->SpL_RE_System    , Box->SmL_On   , 0.5*Box->I_xy*Box->SSD_LLLR, NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLR(Basis_Onsite, A_Basis, Inv, Box->SmL_RE_System    , Box->SpL_On   , 0.5*Box->I_xy*Box->SSD_LLLR, NULL       , &elem_num, NULL   , "No" );
   
   int num,LL,LR,inv;
   for (num = 0; num < elem_num; num++) {
      LL = A_Basis->LL_LLLRRRRL[num];
      LR = A_Basis->LR_LLLRRRRL[num];
      inv = Inv[LL][LR];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
