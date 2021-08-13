//
//  MAKE_ELEMENT_HAM_LLLRRRRL.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
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
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_RE_System  , Box->CUp_On    , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_RE_System    , Box->CUp_D_On  , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_RE_System, Box->CDown_On  , Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LL_LR", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_RE_System  , Box->CDown_D_On, Box->t1*Box->SSD_LLLR, Box->Ele_LL, &elem_num, "LR_LL", "Yes");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_RE_System  , Box->NC_Up_On  , Box->V*Box->SSD_LLLR , NULL       , &elem_num, NULL   , "No" );
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_RE_System, Box->NC_Down_On, Box->V*Box->SSD_LLLR , NULL       , &elem_num, NULL   , "No" );

   //Interaction RRRL
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_RE_Enviro  , Box->CUp_On    , Box->t1*Box->SSD_RRRL, Box->Ele_RR, &elem_num, "RR_RL", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_RE_Enviro    , Box->CUp_D_On  , Box->t1*Box->SSD_RRRL, Box->Ele_RR, &elem_num, "RL_RR", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_RE_Enviro, Box->CDown_On  , Box->t1*Box->SSD_RRRL, Box->Ele_RR, &elem_num, "RR_RL", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_RE_Enviro  , Box->CDown_D_On, Box->t1*Box->SSD_RRRL, Box->Ele_RR, &elem_num, "RL_RR", "Yes");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_RE_Enviro  , Box->NC_Up_On  , Box->V*Box->SSD_RRRL , NULL       , &elem_num, NULL   , "No");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_RE_Enviro, Box->NC_Down_On, Box->V*Box->SSD_RRRL , NULL       , &elem_num, NULL   , "No");
   
   
   if (strcmp(Box->BC, "OBC") == 0 || strcmp(Box->BC, "SSD") == 0 || strcmp(Box->BC, "PBC_LL_LR_RL_RR") == 0) {
      //Interaction LRRL
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_On  , Box->CUp_On    , Box->t1*Box->SSD_LRRL, Box->Ele_RR, Box->Ele_On, &elem_num, "LR_RL", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_On    , Box->CUp_D_On  , Box->t1*Box->SSD_LRRL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LR", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_On, Box->CDown_On  , Box->t1*Box->SSD_LRRL, Box->Ele_RR, Box->Ele_On, &elem_num, "LR_RL", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_On  , Box->CDown_D_On, Box->t1*Box->SSD_LRRL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LR", "Yes");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_On  , Box->NC_Up_On  , Box->V*Box->SSD_LRRL , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_On, Box->NC_Down_On, Box->V*Box->SSD_LRRL , NULL       , NULL       , &elem_num, NULL   , "No" );
   }
   
   if ( strcmp(Box->BC, "PBC_LL_LR_RL_RR") == 0) {
      //Interaction LLRR
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_LE_System  , Box->CUp_LE_Enviro    , Box->t1, Box->Ele_LL, Box->Ele_On, &elem_num, "LL_RR", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_LE_System    , Box->CUp_D_LE_Enviro  , Box->t1, Box->Ele_LL, Box->Ele_On, &elem_num, "RR_LL", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_LE_System, Box->CDown_LE_Enviro  , Box->t1, Box->Ele_LL, Box->Ele_On, &elem_num, "LL_RR", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_LE_System  , Box->CDown_D_LE_Enviro, Box->t1, Box->Ele_LL, Box->Ele_On, &elem_num, "RR_LL", "Yes");
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_LE_System  , Box->NC_Up_LE_Enviro  , Box->V , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_LE_System, Box->NC_Down_LE_Enviro, Box->V , NULL       , NULL       , &elem_num, NULL   , "No" );
   }
   
   if (strcmp(Box->BC, "PBC_LL_LR_RR_RL") == 0) {
      //Interaction LRRR
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_On  , Box->CUp_LE_Enviro    , Box->t1, Box->Ele_On, &elem_num, "LR_RR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_On    , Box->CUp_D_LE_Enviro  , Box->t1, Box->Ele_On, &elem_num, "RR_LR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_On, Box->CDown_LE_Enviro  , Box->t1, Box->Ele_On, &elem_num, "LR_RR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_On  , Box->CDown_D_LE_Enviro, Box->t1, Box->Ele_On, &elem_num, "RR_LR", "Yes");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_On  , Box->NC_Up_LE_Enviro  , Box->V , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_On, Box->NC_Down_LE_Enviro, Box->V , NULL       , &elem_num, NULL   , "No" );

      //Interaction LLRL
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_D_LE_System  , Box->CUp_On    , Box->t1, Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "LL_RL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CUp_LE_System    , Box->CUp_D_On  , Box->t1, Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_D_LE_System, Box->CDown_On  , Box->t1, Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "LL_RL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->CDown_LE_System  , Box->CDown_D_On, Box->t1, Box->Ele_LL, Box->Ele_RR, Box->Ele_On, &elem_num, "RL_LL", "Yes");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Up_LE_System  , Box->NC_Up_On  , Box->V , NULL       , NULL       , NULL       , &elem_num, NULL   , "No" );
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->NC_Down_LE_System, Box->NC_Down_On, Box->V , NULL       , NULL       , NULL       , &elem_num, NULL   , "No" );
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
