//
//  GET_BLOCK_MATRIX.c
//  1DHUBBARD_VF_DMRG
//
//  Created by Kohei Suzuki on 2020/01/24.
//  Copyright © 2020 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

BLOCK *GET_BLOCK_MATRIX(MODEL_1DHUBBARD_VF *Model, DMRG_PARAMETER *Dmrg_Param) {
   
   BLOCK *Block = malloc(sizeof(*Block));
   
   int tot_site   = Model->tot_site;
   int max_dim    = Dmrg_Param->max_dim_system;
   int dim_onsite = Model->dim_onsite;
   int elem_num   = max_dim*max_dim*Dmrg_Param->sp_LL;
   
   Block->Ham          = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_RE       = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_RE     = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CUp_D_RE     = GET_CRS2(tot_site, max_dim, elem_num);
   Block->CDown_D_RE   = GET_CRS2(tot_site, max_dim, elem_num);
   Block->NC_Up_RE     = GET_CRS2(tot_site, max_dim, elem_num);
   Block->NC_Down_RE   = GET_CRS2(tot_site, max_dim, elem_num);

   if (strcmp(Model->BC, "PBC_LL_LR_RR_RL") == 0 || strcmp(Model->BC, "PBC_LL_LR_RL_RR") == 0) {
      Block->CUp_LE       = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_LE     = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CUp_D_LE     = GET_CRS2(tot_site, max_dim, elem_num);
      Block->CDown_D_LE   = GET_CRS2(tot_site, max_dim, elem_num);
      Block->NC_Up_LE     = GET_CRS2(tot_site, max_dim, elem_num);
      Block->NC_Down_LE   = GET_CRS2(tot_site, max_dim, elem_num);
   }
   
   Block->Tot_Sz  = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Tot_Ele = GET_ARRAY_INT2(tot_site, max_dim);
   Block->Dim     = GET_ARRAY_INT1(tot_site);
   
   Block->CUp_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CUp_D_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->CDown_D_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_Up_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_Down_On = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->Ham_On     = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   
   ONSITE_CUP_SZBASIS_HUBBARD           (Block->CUp_On    , 1.0);
   ONSITE_CDOWN_SZBASIS_HUBBARD         (Block->CDown_On  , 1.0);
   ONSITE_CUP_DAGGER_SZBASIS_HUBBARD    (Block->CUp_D_On  , 1.0);
   ONSITE_CDOWN_DAGGER_SZBASIS_HUBBARD  (Block->CDown_D_On, 1.0);
   ONSITE_NC_UP_SZBASIS_HUBBARD         (Block->NC_Up_On  , 1.0);
   ONSITE_NC_DOWN_SZBASIS_HUBBARD       (Block->NC_Down_On, 1.0);
   ONSITE_MAKE_HAM_SZBASIS_HUBBARD      (Model, Block->Ham_On);
   
   //For Expectation Values
   Block->SzC_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->SxC_On  = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->NC_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);
   Block->DO_On   = GET_CRS1(dim_onsite, dim_onsite*dim_onsite);

      
   ONSITE_SZC_SZBASIS_HUBBARD(Block->SzC_On , 1.0);
   ONSITE_SXC_SZBASIS_HUBBARD(Block->SxC_On , 1.0);
   ONSITE_NC_SZBASIS_HUBBARD (Block->NC_On  , 1.0);
   ONSITE_NC_UP_NC_DOWN_SZBASIS_HUBBARD(Block->DO_On  , 1.0);

   Block->TM   = malloc(sizeof(CCS1*)*tot_site/2);
   Block->TM_D = malloc(sizeof(CRS1*)*tot_site/2);
   
   Block->Basis_LL_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_LR_LLLR  = malloc(sizeof(short*)*tot_site/2);
   Block->Basis_Inv_LLLR = GET_ARRAY_INT3(tot_site/2, max_dim, dim_onsite);
   Block->Dim_LLLR       = GET_ARRAY_INT1(tot_site/2);
   
   return Block;
   
}
