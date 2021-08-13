//
//  MAKE_ELEMENT_HAM_LLLRRRRL.c
//  1DKLM_TVF_DMRG
//
//  Created by Kohei Suzuki on 2019/07/15.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box) {
   
   int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_LL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_System, 1.0        , &elem_num);
   DMRG_MAKE_ELEM_RR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_Enviro, 1.0        , &elem_num);
   DMRG_MAKE_ELEM_LR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , Box->SSD_LR, &elem_num);
   DMRG_MAKE_ELEM_RL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , Box->SSD_RL, &elem_num);
   
   //Interaction LLLR
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_RE_System , Box->Even_On   , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_RE_System   , Box->Even_D_On , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_RE_System  , Box->Odd_On    , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_RE_System    , Box->Odd_D_On  , Box->t*Box->SSD_LLLR       , Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_RE_System    , Box->SzL_On    , Box->I_z*Box->SSD_LLLR     , NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_RE_System    , Box->SmL_On    , 0.5*Box->I_xy*Box->SSD_LLLR, NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_RE_System    , Box->SpL_On    , 0.5*Box->I_xy*Box->SSD_LLLR, NULL       , &elem_num, NULL   , "No" );
   
   //Interaction RRRL
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_RE_Enviro , Box->Even_On   , Box->t*Box->SSD_RRRL       , Box->Ele_RR, &elem_num, "RR_RL", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_RE_Enviro   , Box->Even_D_On , Box->t*Box->SSD_RRRL       , Box->Ele_RR, &elem_num, "RL_RR", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_RE_Enviro  , Box->Odd_On    , Box->t*Box->SSD_RRRL       , Box->Ele_RR, &elem_num, "RR_RL", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_RE_Enviro    , Box->Odd_D_On  , Box->t*Box->SSD_RRRL       , Box->Ele_RR, &elem_num, "RL_RR", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_RE_Enviro    , Box->SzL_On    , Box->I_z*Box->SSD_RRRL     , NULL       , &elem_num, NULL   , "No");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_RE_Enviro    , Box->SmL_On    , 0.5*Box->I_xy*Box->SSD_RRRL, NULL       , &elem_num, NULL   , "No");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_RE_Enviro    , Box->SpL_On    , 0.5*Box->I_xy*Box->SSD_RRRL, NULL       , &elem_num, NULL   , "No");
   
   if (strcmp(Box->BC, "OBC") == 0 || strcmp(Box->BC, "SSD") == 0 || strcmp(Box->BC, "PBC_LL_LR_RL_RR") == 0) {
      //Interaction LRRL
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_On, Box->Even_On  , Box->t*Box->SSD_LRRL       , Box->Ele_RR, Box->Ele_On, &elem_num, "LR_RL", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_On  , Box->Even_D_On, Box->t*Box->SSD_LRRL       , Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LR", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_On , Box->Odd_On   , Box->t*Box->SSD_LRRL       , Box->Ele_RR, Box->Ele_On, &elem_num, "LR_RL", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_On   , Box->Odd_D_On , Box->t*Box->SSD_LRRL       , Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LR", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_On   , Box->SzL_On   , Box->I_z*Box->SSD_LRRL     , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_On   , Box->SmL_On   , 0.5*Box->I_xy*Box->SSD_LRRL, NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_On   , Box->SpL_On   , 0.5*Box->I_xy*Box->SSD_LRRL, NULL       , NULL       , &elem_num, NULL   , "No" );
   }
   
   if ( strcmp(Box->BC, "PBC_LL_LR_RL_RR") == 0) {
      //Interaction LLRR
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_LE_System, Box->Even_LE_Enviro  , Box->t       , Box->Ele_LL, Box->Ele_On, &elem_num, "LL_RR", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_LE_System  , Box->Even_D_LE_Enviro, Box->t       , Box->Ele_LL, Box->Ele_On, &elem_num, "RR_LL", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_LE_System , Box->Odd_LE_Enviro   , Box->t       , Box->Ele_LL, Box->Ele_On, &elem_num, "LL_RR", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_LE_System   , Box->Odd_D_LE_Enviro , Box->t       , Box->Ele_LL, Box->Ele_On, &elem_num, "RR_LL", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_LE_System   , Box->SzL_LE_Enviro   , Box->I_z     , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_LE_System   , Box->SmL_LE_Enviro   , 0.5*Box->I_xy, NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_LE_System   , Box->SpL_LE_Enviro   , 0.5*Box->I_xy, NULL       , NULL       , &elem_num, NULL   , "No" );
   }
   
   if (strcmp(Box->BC, "PBC_LL_LR_RR_RL") == 0) {
      //Interaction LRRR
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_On, Box->Even_LE_Enviro  , Box->t       , Box->Ele_On, &elem_num, "LR_RR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_On  , Box->Even_D_LE_Enviro, Box->t       , Box->Ele_On, &elem_num, "RR_LR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_On , Box->Odd_LE_Enviro   , Box->t       , Box->Ele_On, &elem_num, "LR_RR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_On   , Box->Odd_D_LE_Enviro , Box->t       , Box->Ele_On, &elem_num, "RR_LR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_On   , Box->SzL_LE_Enviro   , Box->I_z     , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_On   , Box->SmL_LE_Enviro   , 0.5*Box->I_xy, NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_On   , Box->SpL_LE_Enviro   , 0.5*Box->I_xy, NULL       , &elem_num, NULL   , "No" );
      
      //Interaction LLRL
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_D_LE_System, Box->Even_On  , Box->t       , Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "LL_RL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Even_LE_System  , Box->Even_D_On, Box->t       , Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_D_LE_System , Box->Odd_On   , Box->t       , Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "LL_RL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Odd_LE_System   , Box->Odd_D_On , Box->t       , Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SzL_LE_System   , Box->SzL_On   , Box->I_z     , NULL       , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SpL_LE_System   , Box->SmL_On   , 0.5*Box->I_xy, NULL       , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->SmL_LE_System   , Box->SpL_On   , 0.5*Box->I_xy, NULL       , NULL       , NULL       , &elem_num, NULL   , "No" );
   }
   
   DMRG_MAKE_ELEM_ZERO_LLLRRRRL(Basis_Onsite, A_Basis, Inv, &elem_num);
   
   int dim_RR     = Basis_Onsite->dim_RR;
   int dim_onsite = Basis_Onsite->dim_onsite;
   int num,LL,LR,RR,RL,inv;
   long inv_sup;
   for (num = 0; num < elem_num; num++) {
      LL      = A_Basis->LL_LLLRRRRL[num];
      LR      = A_Basis->LR_LLLRRRRL[num];
      RL      = A_Basis->RL_LLLRRRRL[num];
      RR      = A_Basis->RR_LLLRRRRL[num];
      inv_sup = (long)LL*dim_onsite*dim_RR*dim_onsite + (long)LR*dim_RR*dim_onsite + RR*dim_onsite + RL;
      inv     = Inv[inv_sup];
      A_Basis->Inv_LLLRRRRL[inv] = -1;
   }
   
   A_Basis->elem_num = elem_num;
   
}
