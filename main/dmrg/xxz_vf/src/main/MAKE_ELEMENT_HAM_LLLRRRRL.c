//
//  MAKE_ELEMENT_HAM_LLLRRRRL.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/14.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

void MAKE_ELEMENT_HAM_LLLRRRRL(DMRG_BASIS_ONSITE *Basis_Onsite, DMRG_A_BASIS *A_Basis, int *Inv, HAM_BOX *Box) {
   
   int elem_num = 0;
   
   //Onsite Ham
   DMRG_MAKE_ELEM_LL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_System, 1.0, &elem_num);
   DMRG_MAKE_ELEM_RR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_Enviro, 1.0, &elem_num);
   DMRG_MAKE_ELEM_LR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0, &elem_num);
   DMRG_MAKE_ELEM_RL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Ham_On    , 1.0, &elem_num);
   
   //Interaction LLLR
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_RE_System, Box->Sm_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_RE_System, Box->Sp_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_LLLR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_RE_System, Box->Sz_On, Box->J_z , NULL, &elem_num, NULL, "No");
   
   //Interaction RRRL
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_RE_Enviro, Box->Sm_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_RE_Enviro, Box->Sp_On, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
   DMRG_MAKE_ELEM_RRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_RE_Enviro, Box->Sz_On, Box->J_z , NULL, &elem_num, NULL, "No");
   
   if (strcmp(Box->BC, "OBC") == 0 || strcmp(Box->BC, "SSD") == 0) {
      //Interaction LRRL
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_On, Box->Sm_On, 0.5*Box->J_xy, NULL, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_On, Box->Sp_On, 0.5*Box->J_xy, NULL, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LRRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_On, Box->Sz_On, Box->J_z , NULL, NULL, &elem_num, NULL, "No");
   }
   
   else if (strcmp(Box->BC, "PBC") == 0) {
      //Interaction LRRR
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_On, Box->Sm_LE_Enviro, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_On, Box->Sp_LE_Enviro, 0.5*Box->J_xy, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LRRR_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_On, Box->Sz_LE_Enviro, Box->J_z , NULL, &elem_num, NULL, "No");
      
      //Interaction LLRL
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sp_LE_System, Box->Sm_On, 0.5*Box->J_xy, NULL, NULL, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sm_LE_System, Box->Sp_On, 0.5*Box->J_xy, NULL, NULL, NULL, &elem_num, NULL, "No");
      DMRG_MAKE_ELEM_LLRL_LLLRRRRL(Basis_Onsite, A_Basis, Inv, Box->Sz_LE_System, Box->Sz_On, Box->J_z , NULL, NULL, NULL, &elem_num, NULL, "No");
   }
   else {
      printf("Error in MAKE_ELEMENT_HAM_LLLRRRRL\n");
      exit(1);
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
